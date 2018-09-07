/*
 * Copyright (c) 2018
 * inviwo-convert is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * inviwo-convert is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with inviwo-convert; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 */
/*! \internal \file
 * \brief
 *
 * \author Christian Blau <cblau@gwdg.de>
 */
#include "mrcfile.h"
#include "mrcheader.h"

#include <cstdio>

#include <algorithm>
#include <complex>
#include <string>
#include <vector>
#include <set>
#include <type_traits>

/*******************************************************************************
 * MrcFileView::Impl
 */
class MrcFileView::Impl
{
    public:
        Impl();
        ~Impl();

        void read_file_size();

        void read_mrc_data_();
        void read_mrc_header_();

        /*! \brief Guess, whether endianess differs between input file and reading architecture .
         *
         * If the number of columns in the density file is negative or larger than 65534,
         * assume endianess missmatch between input file and reading machine architecture.*/
        void check_swap_bytes();

        bool has_skew_matrix();

        template <typename T> void read(T * result)
        {
            fread(result, sizeof(T), 1, file_);
            // swap bytes for correct endianness
            if (header_.swap_bytes)
            {
                // byte swap for real numbers
                if (std::is_same<T, double>())
                {
                    int64_t int_tmp = int64_t(*result);
                    *result = (int_tmp & 0xFF00000000000000) >> 7*8 | (int_tmp & 0x00FF000000000000) >> 5*8 |
                        (int_tmp & 0x0000FF0000000000) >> 3*8 | (int_tmp & 0x000000FF00000000) >> 1*8 |
                        (int_tmp & 0x00000000FF000000) << 1*8 | (int_tmp & 0x0000000000FF0000) << 3*8 |
                        (int_tmp & 0x000000000000FF00) << 5*8 | (int_tmp & 0x00000000000000FF) << 7*8;
                }

                // byte swap for real numbers
                if (std::is_same<T, float>())
                {
                    int32_t int_tmp = int32_t(*result);
                    *result = (int_tmp & 0xFF000000) >> 24 | (int_tmp & 0x00FF0000) >> 8 | (int_tmp & 0x0000FF00) << 8 | (int_tmp & 0x000000FF) << 24;
                }
                if (std::is_same<T, int16_t>())
                {
                    int16_t int_tmp = int16_t(*result);
                    *result = (int_tmp & 0xFF00) >> 8 | (int_tmp & 0x00FF) << 8;
                }
            }

        }

        void read_float32_rvec_(std::array<float,3> * result);
        void read_int32_ivec_(std::array<int, 3> * result);

        bool colummn_row_section_order_valid_(std::array<int, 3> crs_to_xyz);

        FILE                          *file_;
        size_t                         file_size_;
        constexpr static size_t        numLabels_c   = 10;
        constexpr static size_t        labelSize_c   = 80;
        constexpr static size_t        headerBytes_c = 1024;
        constexpr static size_t XX = 0;
        constexpr static size_t YY = 1;
        constexpr static size_t ZZ = 2;


        const std::vector<std::string> filetypes;

        MrcHeader           header_;
        std::vector<float>  data_;

};


void MrcFileView::Impl::read_int32_ivec_(std::array<int, 3> * result)
{
    read(&(*result)[XX]);
    read(&(*result)[YY]);
    read(&(*result)[ZZ]);

}

void MrcFileView::Impl::read_float32_rvec_(std::array<float,3> * result)
{
    read(&(*result)[XX]);
    read(&(*result)[YY]);
    read(&(*result)[ZZ]);
}


bool MrcFileView::Impl::colummn_row_section_order_valid_(std::array<int, 3> crs_to_xyz)
{
    const std::set<int> valid_crs_set {
        0, 1, 2
    };
    std::set<int> crs_set {
        crs_to_xyz[XX], crs_to_xyz[YY], crs_to_xyz[ZZ]
    };
    return valid_crs_set == crs_set;
};

void MrcFileView::Impl::read_mrc_header_()
{
        check_swap_bytes();
        read_file_size();

    /* Supports reading according to
       ftp://ftp.wwpdb.org/pub/emdb/doc/Map-format/current/EMDB_map_format.pdf
       note, that
       http://www.ccpem.ac.uk/mrc_format/mrc2014.php
       differs slightly in definition */

    /* 1-3 | NC, NR, NS | signed int >0
     * # of columns (fastest changing),rows, sections (slowest changing)
     * emdb convention: NC=NR=NS                     */

    read_int32_ivec_(&header_.num_crs);

    /* 4   | MODE | signed int | 0,1,2,3,4
     * voxel datatype
     * emdb convention: 2       */
    read(&header_.mrc_data_mode);

    /* MODE = 0: 8 bits, density stored as a signed byte (range -128 to 127, ISO/IEC 10967)
     * MODE = 1: 16 bits, density stored as a signed integer (range -32768 to 32767, ISO/IEC 10967)
     * MODE = 2: 32 bits, density stored as a floating point number (IEEE 754)
     * MODE = 3: 32 bits, Fourier transform stored as complex signed integers (ISO/IEC 10967)
     * MODE = 4: 64 bits, Fourier transform stored as complex floating point numbers (IEEE 754)     */

    /* 5-7 | NCSTART, NRSTART, NSSTART | signed int
     * position of first column, first row, and first section (voxel grid units)
     *
     * The position of the first voxel is defined in grid units by NCSTART, NRSTART, and NSSTART.
     * The center of the voxel with grid position (0,0,0) corresponds to the Cartesian coordinate origin.*/
    read_int32_ivec_(&header_.crs_start);

    /* 8-10 | NX, NY, NZ | signed int >0 |
     * intervals per unit cell repeat along X,Y Z
     * intervals per map length along X,Y,Z;
     * emdb convention: same as NC, NR, NS
     *
     * Lengths in Aangstroms for a single voxel are as follows:
     * Xvoxel = X_LENGTH/NX Yvoxel = Y_LENGTH/NY Zvoxel = Z_LENGTH/NZ */

    read_int32_ivec_(&header_.extend);

    /* 11-13 | X_LENGTH, Y_LENGTH, Z_LENGTH | floating pt >0
     * Unit Cell repeats along X, Y, Z In Aangstrom
     * emdb Map lengths along X,Y,Z in Aangstrom
     *
     * Lengths in Aangstroms for a single voxel are as follows:
     * Xvoxel = X_LENGTH/NX Yvoxel = Y_LENGTH/NY Zvoxel = Z_LENGTH/NZ */
    read_float32_rvec_(&header_.cell_length);

    /* 14-16 | ALPHA,BETA,GAMMA | floating pt >0, <180
     * Unit Cell angles (degrees)
     * emdb convention: 90, 90, 90
     *
     * By convention, cell angles (ALPHA, BETA, GAMMA)
     * are 90 degrees for single particle or tomogram EM maps;
     * they follow IUCr space group conventions for crystals.*/
    read_float32_rvec_(&header_.cell_angles);
    // By convention, unset cell angles (all 0) are interpreted as 90 deg.
    if (header_.cell_angles[XX]*header_.cell_angles[YY]*header_.cell_angles[ZZ] < 1e-5)
    {
        header_.cell_angles = {90, 90, 90};
    }

    /* 17-19 | MAPC, MAPR, MAPS | signed int | 1 (=X) 2 (=Y) 3 (=Z)
     * relationship of X,Y,Z axes to columns, rows, sections
     * emdb convention: 1, 2, 3 */
    std::array<int, 3> crs_to_xyz {{
                                       header_.crs_to_xyz[XX]+1, header_.crs_to_xyz[YY]+1, header_.crs_to_xyz[ZZ]+1
                                   }};
    read_int32_ivec_(&crs_to_xyz);

        header_.crs_to_xyz[XX] = crs_to_xyz[XX]-1;
        header_.crs_to_xyz[YY] = crs_to_xyz[YY]-1;
        header_.crs_to_xyz[ZZ] = crs_to_xyz[ZZ]-1;
        if (!colummn_row_section_order_valid_(header_.crs_to_xyz))
        {
            header_.crs_to_xyz = {{0, 1, 2}};
        }


    /* 20-22 | AMIN, AMAX, AMEAN | floating pt
     * Minimum, maximum, average density */
    read(&(header_.min_value ));
    read(&(header_.max_value ));
    read(&(header_.mean_value));

    /* 23 | ISPG | signed int 1-230 |
     * space group #
     * emdb convention 1
     *
     * Space Group Numbers are defined by IUCr conventions
     * (Table 12.3.4.1 Standard space-group symbols”, pages 824-831,
     * International Tables for Crystallography, Volume A, fifth edition).
     *
     * For 3D volumes of single particle or tomogram entries, ISPG=1 and NSYMBT=0.
     * For image stacks ISPG = 0 */
    read(&header_.space_group);

    /* 24 | NSYMBT | signed int | 80n
     * # of bytes in symmetry table (multiple of 80)
     * emdb convention 0 */
    int numBytesExtendedHeader;
    read(&numBytesExtendedHeader);

    if (header_.is_crystallographic)
    {
        /* 25 | LSKFLG | signed int | 0,1
         * flag for skew matrix
         * emdb convention 0 */
        int32_t hasSkewMatrix = header_.has_skew_matrix ? 1 : 0;
        read(&hasSkewMatrix);
        header_.has_skew_matrix = (hasSkewMatrix == 1) ? true : false;

        if (header_.has_skew_matrix)
        {
            /* TODO: A2NM conversion for skew matrix if necessary */
            /* 26-34 | SKWMAT | floating pt
             * skew matrix-S11, S12, S13, S21, S22, S23, S31, S32, S33
             * emdb convention: not set
             *
             * SKWMAT, SKWTRN, and EXTRA fields are not currently used by EMDB. */

            /* 35-37 | SKWTRN | floating pt
             * skew translation-T1, T2, T3
             * emdb convention: not set
             *
             * SKWMAT, SKWTRN, and EXTRA fields are not currently used by EMDB. */
            for (auto && i : header_.skew_matrix)
            {
                read(&i);
            }
            read_float32_rvec_(&header_.skew_translation);
        }
    }
    else
    {
        /* 25-37 not used in EMDB */
        for (auto && i : header_.extraskew)
        {
            read(&i);
        }
    }

    /* 38-52 | EXTRA | 32 bit binary
     * user-defined metadata
     *
     * SKWMAT, SKWTRN, and EXTRA fields are not currently used by EMDB.
     * EMDB might use fields 50,51 and 52 for setting the coordinate system origin */
    for (auto && i : header_.extra)
    {
        read(&i);
    }

    /* 53 | MAP | ASCII char
     * "MAP "
     * MRC/CCP4 MAP format identifier */
    header_.format_identifier.resize(4);
    for (size_t i = 0; i < 4; i++)
    {
        read(&header_.format_identifier[i]);
    }
    /* 54 | MACHST | 32 bit
     * binary machine stamp
     *
     * MACHST is (written/read as 4 hex byte sequence)
     * 0x44,0x41,0x00,0x00  for little endian machines
     * 0x11,0x11,0x00,0x00  for big endian machines
     */
    read(&(header_.machine_stamp));

    /* 55 | RMS | floating pt
     * Density root-mean-square deviation */
    read(&(header_.rms_value));

    /* 56 | NLABL | signed int | 0-10
     * # of labels
     *
     * Following the 2010 remediation, maps distributed by EMDB
     * now have a single label of form “::::EMDataBank.org::::EMD-1234::::”.  */
    read(&header_.num_labels);

    /* 57-256 | LABEL_N | ASCII char
     * 10 user-defined labels each 80 characters long */

        for (auto &label : header_.labels)
        {
            std::array<char, labelSize_c>  labelToRead;
            fread(labelToRead.data(), 1, labelSize_c, file_);
            label = std::string(labelToRead.data(), labelSize_c);
        }

    /* 257-257+NSYMBT | anything
     */
    header_.extended_header.resize(numBytesExtendedHeader);

};

void MrcFileView::Impl::read_mrc_data_()
{
    data_.resize(header_.extend[XX] * header_.extend[YY] * header_.extend[ZZ]);

    for (float & d : data_)
    {
        read(&d);
    }
}

void MrcFileView::Impl::read_file_size()
{
    fpos_t current;
    fgetpos(file_, &current);
    fseek(file_, 0, SEEK_END);
    file_size_ = ftell(file_);
    fsetpos(file_, &current);
}

void MrcFileView::Impl::check_swap_bytes()
{

    fpos_t      current;
    header_.swap_bytes = false;
    int32_t number_columns;
    fgetpos(file_, &current);

    fseek(file_, 0, SEEK_SET);
    read(&number_columns);
    if (number_columns <= 0 || number_columns >= 65536)
    {
        header_.swap_bytes = true;
    }

    // rewind the file
    fsetpos(file_, &current);

}

MrcFileView::Impl::Impl() : file_(nullptr), file_size_(0), filetypes({"mrc", "ccp4", "imod", "map"}
                                                                 )
{
    header_.setEMDBDefaults();
};

MrcFileView::Impl::~Impl()
{
    if (file_ != nullptr)
    {
        fclose(file_);
    }
};



/*******************************************************************************
 * MrcFileView
 */

MrcFileView::MrcFileView(const std::string & filename):
impl_(new MrcFileView::Impl)
{
    impl_->file_ = fopen(filename.c_str(), "r");
    impl_->read_mrc_header_();
    impl_->read_mrc_data_();
}

MrcFileView::~MrcFileView()
{

}
const MrcHeader & MrcFileView::header() const
{
    return impl_->header_;
}

const std::vector<float> & MrcFileView::data() const
{
    return impl_->data_;
}
