/**
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// TCLAP
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "BaseLib/LogogSimpleFormatter.h"

// FileIO
#include "FileIO/readMeshFromFile.h"
#include "FileIO/XmlIO/Boost/BoostVtuInterface.h"

// MathLib
#include "MathLib/Vector3.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshSurfaceExtraction.h"

int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format(new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Extract surface of given mesh.", ' ', "0.1");
	TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
		"the name of the file containing the input mesh",
		true,
		"",
		"file name as string");
	cmd.add(mesh_in);
	TCLAP::ValueArg<std::string> mesh_out("o", "mesh-output-file",
		"the name of the file the surface mesh will be written to",
		true,
		"",
		"file name as string");
	cmd.add(mesh_out);
	cmd.parse(argc, argv);

	MeshLib::Mesh mesh(*FileIO::readMeshFromFile(mesh_in.getValue()));
	INFO("Mesh read: %d nodes, %d elements.", mesh.getNNodes(), mesh.getNElements());

	 // this extracts the complete surface
	MathLib::Vector3 const dir(0.0,0.0,0.0);
	MeshLib::Mesh surface_mesh(
		*MeshLib::MeshSurfaceExtraction::getMeshSurface(
			mesh,
			dir,
			90)
		);

	// write surface mesh into a file
	FileIO::BoostVtuInterface out;
	out.setMesh(&surface_mesh);
	out.writeToFile(mesh_out.getValue());

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}

