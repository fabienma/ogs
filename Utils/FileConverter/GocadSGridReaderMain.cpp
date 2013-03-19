/**
 * This file is part of "GocadSGridReader". "GocadSGridReader" is free
 * software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * "GocadSGridReader" is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License with
 * "GocadSGridReader". If not, see <http://www.gnu.org/licenses/>.
 *
 * @date 2013-03-01
 * @author Thomas Fischer
 * @file main.cpp
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "tclap/CmdLine.h"
#include "FileTools.h"
#include "LogogSimpleFormatter.h"

// FileIO
#include "RapidXmlIO/BoostVtuInterface.h"

// MeshLib
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"

// Utils/FileConverter
#include "GocadSGridReader.h"

int main(int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog::Cout *logog(new logog::Cout);
	logog->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Programme reads parts of a Gocad structured mesh", ' ', "0.1");

	// Define a value argument and add it to the command line.
	// A value arg defines a flag and a type of value that it expects
	TCLAP::ValueArg<std::string> sg_file_arg("s", "sg", "structured grid file name", true, "",
			"string");

	// Add the argument sg_file_arg to the CmdLine object. The CmdLine object
	// uses this Arg to parse the command line.
	cmd.add(sg_file_arg);

	cmd.parse(argc, argv);

	// read the Gocad SGrid
	FileIO::GocadSGridReader reader(sg_file_arg.getValue());
	std::vector<MeshLib::Node*> const& nodes(reader.getNodes());
	std::vector<MeshLib::Element*> const& elements(reader.getElements());

	MeshLib::Mesh mesh("GocadSGrid", nodes, elements);

	FileIO::BoostVtuInterface vtu;
	vtu.setMesh(&mesh);
	// output file name
	std::string mesh_out_fname(BaseLib::extractPath(sg_file_arg.getValue()) +
				BaseLib::extractBaseNameWithoutExtension(sg_file_arg.getValue()) + ".vtu");
	vtu.writeToFile(mesh_out_fname);

	delete logog;
	delete custom_format;
	LOGOG_SHUTDOWN();

	return 0;
}
