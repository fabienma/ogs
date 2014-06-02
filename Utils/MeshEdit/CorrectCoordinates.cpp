/**
 * @file CorrectCoordinates.cpp
 * @date 2014-06-02
 * @brief
 *
 * @copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// ThirdParty
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "StringTools.h"
#include "FileTools.h"
#include "LogogSimpleFormatter.h"

// FileIO
#include "readMeshFromFile.h"
#include "XmlIO/Boost/BoostVtuInterface.h"

// GeoLib
#include "AABB.h"

// MeshLib
#include "Node.h"
#include "Elements/Element.h"
#include "MeshEditing/moveMeshNodes.h"
#include "Mesh.h"

int main(int argc, char *argv[])
{
	LOGOG_INITIALIZE();
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog::Cout *logogCout(new logog::Cout);
	logogCout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Moves the mesh nodes using the given displacement vector or if no displacement vector is given, moves the mesh nodes such that the centroid of the given mesh is in the origin.", ' ', "0.1");
	// Define a value argument and add it to the command line.
	// A value arg defines a flag and a type of value that it expects,
	// such as "-m meshfile".
	TCLAP::ValueArg<std::string> mesh_arg("m","mesh","input mesh file",true,"","string");

	// Add the argument mesh_arg to the CmdLine object. The CmdLine object
	// uses this Arg to parse the command line.
	cmd.add( mesh_arg );

	TCLAP::ValueArg<std::string> mesh_out_arg("o","output-mesh","output mesh file", false, "", "string");
	cmd.add(mesh_out_arg);

	cmd.parse( argc, argv );

	std::string fname (mesh_arg.getValue());

	MeshLib::Mesh* mesh = FileIO::readMeshFromFile(fname);

	std::for_each(mesh->getNodes().begin(), mesh->getNodes().end(),
		[](MeshLib::Node* node) {
			MeshLib::Node &n(*node);
			if (std::abs(n[0]) <= 1e-4)
				n[0] = 0.0;
			if (std::abs(n[1]) <= 1e-4)
				n[1] = 0.0;
			if (std::abs(n[2]) <= 1e-4)
				n[2] = 0.0;
		}
	);

	std::string out_fname(mesh_out_arg.getValue());
	if (out_fname.empty()) {
		out_fname = BaseLib::dropFileExtension(mesh_out_arg.getValue());
		out_fname += "_corrected.vtu";
	}

	FileIO::BoostVtuInterface mesh_io;
	mesh_io.setMesh(mesh);
	mesh_io.writeToFile(out_fname);

	delete mesh;
	delete logogCout;
	delete custom_format;
	LOGOG_SHUTDOWN();
}


