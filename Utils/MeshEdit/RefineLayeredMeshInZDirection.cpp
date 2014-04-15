/**
 * @file RefineLayeredMeshInZDirection.cpp
 * @date Jan 17, 2014
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
#include "Mesh.h"
#include "MeshEnums.h"

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

	MeshLib::Mesh* original_mesh = FileIO::readMeshFromFile(fname);

	// check if mesh does not contain pyramids
	std::vector<MeshLib::Element*> const& original_elements(original_mesh->getElements());
	std::size_t const n_elements(original_elements.size());
	for (std::size_t k(0); k<n_elements; k++) {
		MeshElemType const type(original_elements[k]->getGeomType());
		if (type != MeshElemType::TETRAHEDRON
			&& type != MeshElemType::PRISM
			&& type != MeshElemType::HEXAHEDRON) {
			INFO("Element %d of mesh is not a tet, prism or hex.", k);
			delete original_mesh;
			break;
		}
	}

	if (!original_mesh) {
		return -1;
	}

	std::string out_fname(mesh_out_arg.getValue());
	if (out_fname.empty()) {
		out_fname = BaseLib::dropFileExtension(mesh_out_arg.getValue());
		out_fname += "_refined.vtu";
	}

	FileIO::BoostVtuInterface mesh_io;
	mesh_io.setMesh(original_mesh);
	mesh_io.writeToFile(out_fname);

	delete original_mesh;
	delete logogCout;
	delete custom_format;
	LOGOG_SHUTDOWN();
}


