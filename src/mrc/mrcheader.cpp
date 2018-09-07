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
 * Implements mrc header.
 *
 * \author Christian Blau <cblau@gwdg.de>
 */

#include "mrcheader.h"

const float MrcHeader::nmToMrcUnits = 10.0f;

void MrcHeader::setEMDBDefaults()
{
    {
        swap_bytes               = false;
        space_group              = 1;
        mrc_data_mode            = 2;
        num_bytes_extened_header = 0;
        has_skew_matrix          = false;
        crs_start                = {{0, 0, 0}};
        crs_to_xyz               = {{0, 1, 2}};
        skew_matrix              = {{0, 0, 0, 0, 0, 0, 0, 0, 0}};
        skew_translation         = {0, 0, 0};
        is_crystallographic      = false;
        extra                    = {};
        extraskew                = {};
        format_identifier        = "MAP ";

        machine_stamp            = 1145110528;
        // machine_stamp            = 4369;
        std::string empty80CharLabel = {"                                                                                "};
        std::string emdbCustomLabel  = {"::::EMDataBank.org::::EMD-xxxx::::Own Data Following EMDB convention::::::::::::"};
        labels[0] = emdbCustomLabel;
        for (int i = 1; i < 10; ++i)
        {
            labels[i] = empty80CharLabel;
        }
        num_labels               = 1;
        extended_header          = {};
    }
}
