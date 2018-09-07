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
/*! \file
 * \brief
 * Reading and writing routines for volume data formats ccp4, mrc and imod.
 *
 * \author Christian Blau <cblau@gwdg.de>
 */

#ifndef MRCFILE_H_
#define MRCFILE_H_

#include <memory>
#include <vector>
struct MrcHeader;

 /*! \brief View an Mrc File.
 *
 * Read float valued volume data files
 * according to the electron microscopy data bank (EMDB) standard.
 *
 * The formatting guraranties compliance with 3D EM maps described in
 * "EMDB Map Distribution Format Description Version 1.01 (c) emdatabank.org 2014"
 *
 * However, other ccp4, mrc, imod and map formats might be compatible.
 * \param[in] filename name of the file from which to read the griddata, typically *.cpp4, *.mrc or *.map
 * \returns MrcFileView into float-valued, real-space data on a grid.
 */
class MrcFileView
{
public:
    explicit MrcFileView(const std::string & filename);
    ~MrcFileView();
    const MrcHeader & header() const;
    const std::vector<float> & data() const;
private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

#endif /* end of include guard: MRCFILE_H_ */
