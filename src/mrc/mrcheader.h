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
 *
 */
/*! \file
 * \brief
 * Data structure to hold the complete mrc-file metadata as specified in the mrc format.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 * \ingroup module_fileio
 */

#ifndef MRCHEADER_H_
#define MRCHEADER_H_

#include <array>
#include <string>
#include <vector>
#include <memory>


/*! \brief
 * A container for the metadata in mrc file formats (compatible with ccp4 and map and mostly imod).
 *
 * For a detailed decription see
 * "EMDB Map Distribution Format Description Version 1.01 (c) emdatabank.org 2014"
 */
struct MrcHeader{
    bool                       swap_bytes;               //!< swap bytes upon reading/writing (applied, when endianess is different between file and machine architecture)
    int                        space_group;              //!< space group as defined by IUCr conventions (Table 12.3.4.1 Standard space-group symbols, pages 824-831, International Tables for Crystallography, Volume A, fifth edition)
    /*!\brief The mrc standard defines modes 0-4.
     *
     * MODE = 0: 8 bits, density stored as a signed byte (range -128 to 127, ISO/IEC 10967)
     * MODE = 1: 16 bits, density stored as a signed integer (range -32768 to 32767, ISO/IEC 10967)
     * MODE = 2: 32 bits, density stored as a floating point number (IEEE 754)
     * MODE = 3: 32 bits, Fourier transform stored as complex signed integers (ISO/IEC 10967)
     * MODE = 4: 64 bits, Fourier transform stored as complex floating point numbers (IEEE 754)
     */
    enum class                    MrcDataMode : int { uInt8 = 0, int16 = 1, float32 = 2, complexInt32 = 3, complexFloat64 = 4 };
    int                           mrc_data_mode;            //!< data mode, currently only mode 2 is supported (32-bit float real values)
    int                           machine_stamp;            //!< endianess of map writing architecture (big endian: 0x44410000 , little endian: 0x11110000)
    std::string                   format_identifier;        //!< for all density formats: four 1-byte chars reading "MAP " (a little pointless, I know)

    int                           num_labels;               //!< number of used crystallographic labels, 0 for imagestacks, 1 for emdb data
    std::array<std::string, 10>   labels;                   //!< crystallographic labels or ::::EMDataBank.org::::EMD-1234:::: for EMDB entries

    std::array<float,3>                          cell_length;              //!< length of the crystallographic unit cell
    std::array<float,3>                          cell_angles;              //!< crystallographic unit cell angles

    std::array<int, 3>            crs_to_xyz;               //!< Axis order
    std::array<int, 3>            num_crs;                  //!< redundand entry, we use the grid extend (NX,NY,NZ) from header words 8-10
    std::array<int, 3>            extend;                   //!< the grid extend, check against num_crs
    std::array<int, 3>            crs_start;                //!< Start of values in grid, typically 0,0,0

    float                         min_value;                //!< minimum voxel value. may be used to scale values in currently unsupported compressed data mode (mrc_data_mode=0)
    float                         max_value;                //!< maximum voxel value. may be used to scale values in currently unsupported compressed data mode (mrc_data_mode=0)
    float                         mean_value;               //!< mean voxel value   (not always reported,as evident from density)
    float                         rms_value;                //!< rms of the density (not always reported,as evident from density)

    bool                          is_crystallographic;      //!< true if crystallographic data is to be read
    bool                          has_skew_matrix;          //!< only crystallographic data: true if skew matrix is stored
    std::array<float, 9>          skew_matrix;              //!< only crystallographic data: skew matrix or, if skew flag is zero, data in place of skew matrix
    std::array<float,3>                          skew_translation;         //!< only crystallographic data: skew translatation or, if skew flag is zero, data in place of skew translation
    int                           num_bytes_extened_header; //!< only crystallographic data: the size of the symbol table in bytes
    std::vector<char>             extended_header;          //!< only crystallographic data: extended header, usually symbol tables

    std::array<float, 13>         extraskew;                //!< fields unused in EMDB standard, but used for skew matrix and translation in crystallogrphic data (skew flag, skew matrix and skew translation)
    std::array<float, 15>         extra;                    //!< extra data in header, currently unused

    static const float            nmToMrcUnits;             //!< Conversion factor from nm to mrc units (Angstrom)
    /*! \brief
     * Set values to emdb defauts.
     */
    void setEMDBDefaults();
};

#endif /* end of include guard: MRCHEADER_H_ */
