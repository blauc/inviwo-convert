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
#include <string>
#include <fstream>

#include "mrc/mrcfile.h"
#include "mrc/mrcheader.h"

int main(int argc, const char *argv[]) {

	std::string filename(argv[1]);
	const MrcFileView mrcfile(filename);

	std::string rawFileName = filename + ".raw";
	std::ofstream rawDataStream(rawFileName,  std::ios::out | std::ios::binary);
    rawDataStream.write(reinterpret_cast<const char*>(mrcfile.data().data()), sizeof(float) * mrcfile.data().size());
	rawDataStream.close();
    fprintf(stderr,"Dumped voxel data into \"%s\"\n",rawFileName.c_str());


	const auto & header = mrcfile.header();
	std::string headerFileName = filename + ".dat";
	std::ofstream headerStream(headerFileName,  std::ios::out);
	headerStream << "Rawfile: "  << rawFileName << std::endl;
	headerStream << "Resolution: " << header.extend[0] << " " << header.extend[1] << " " << header.extend[2] << std::endl;
	headerStream << "Format: " << "FLOAT32" << std::endl;
	headerStream << "BasisVector1: " << header.cell_length[0] << " 0 0" << std::endl;
	headerStream << "BasisVector2: " << "0 " << header.cell_length[1] << " 0" << std::endl;
	headerStream << "BasisVector3: " << "0 0 " << header.cell_length[2] << std::endl;
	headerStream.close();
	fprintf(stderr,"Converted header to \"%s\"\n",headerFileName.c_str());

	fprintf(stderr,"Done\n");
	return 0;
}
