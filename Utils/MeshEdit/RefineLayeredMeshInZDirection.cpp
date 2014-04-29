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

#include <algorithm>

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

// MathLib
#include "MathTools.h"

// GeoLib
#include "AABB.h"

// MeshLib
#include "Node.h"
#include "Elements/Element.h"
#include "Elements/Prism.h"
#include "Elements/Tet.h"
#include "Mesh.h"
#include "MeshEnums.h"
#include "MeshEditing/DuplicateMeshComponents.h"

/** Split a prism into two prisms */
void splitElement(MeshLib::Prism* original_prism,
	std::vector<MeshLib::Node*> & nodes,
	std::vector<MeshLib::Element*> & elements)
{
	// deep copy of nodes of original prism
	for (std::size_t k(0); k<6; k++)
		nodes.push_back(new MeshLib::Node(*original_prism->getNodes()[k]));

	std::size_t const s(nodes.size());
	// nodes used for split
	nodes.push_back(new MeshLib::Node(
		0.5 * ((*nodes[s-6])[0] + (*nodes[s-3])[0]),
		0.5 * ((*nodes[s-6])[1] + (*nodes[s-3])[1]),
		0.5 * ((*nodes[s-6])[2] + (*nodes[s-3])[2])
	));
	nodes.push_back(new MeshLib::Node(
		0.5 * ((*nodes[s-5])[0] + (*nodes[s-2])[0]),
		0.5 * ((*nodes[s-5])[1] + (*nodes[s-2])[1]),
		0.5 * ((*nodes[s-5])[2] + (*nodes[s-2])[2])
	));
	nodes.push_back(new MeshLib::Node(
		0.5 * ((*nodes[s-4])[0] + (*nodes[s-1])[0]),
		0.5 * ((*nodes[s-4])[1] + (*nodes[s-1])[1]),
		0.5 * ((*nodes[s-4])[2] + (*nodes[s-1])[2])
	));

	std::array<MeshLib::Node*, 6> nodes_prism0 =
		{{nodes[s-6], nodes[s-5], nodes[s-4], nodes[s], nodes[s+1], nodes[s+2]}};
	std::array<MeshLib::Node*, 6> nodes_prism1 =
		{{nodes[s], nodes[s+1], nodes[s+2], nodes[s-3], nodes[s-2], nodes[s-1]}};
	elements.push_back(new MeshLib::Prism(nodes_prism0, original_prism->getValue()));
	elements.push_back(new MeshLib::Prism(nodes_prism1, original_prism->getValue()));
}

MeshLib::Mesh refineMesh(MeshLib::Mesh const& original_mesh)
{
	std::vector<MeshLib::Node*> nodes;
	std::vector<MeshLib::Element*> elements;

	std::vector<MeshLib::Element*> const& original_elements(original_mesh.getElements());
	std::size_t const n_elements(original_elements.size());
	for (std::size_t k(0); k<n_elements; k++) {
		MeshElemType const type(original_elements[k]->getGeomType());
		switch (type) {
		case MeshElemType::TETRAHEDRON:
			break;
		case MeshElemType::PRISM:
			splitElement(dynamic_cast<MeshLib::Prism*>(original_elements[k]), nodes, elements);
			break;
		case MeshElemType::HEXAHEDRON:
		default: {
			INFO("Element %d of mesh is not a tet, prism or hex.", k);
		}
		}
	}

	return MeshLib::Mesh("RefinedMesh", nodes, elements);
}

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

	double tets_vol(0.0);
	double prisms_vol(0.0);

	for (std::size_t k(0); k<n_elements; k++) {
		MeshElemType const type(original_elements[k]->getGeomType());
		if (type == MeshElemType::TETRAHEDRON)
			tets_vol += original_elements[k]->getContent();
		if (type == MeshElemType::PRISM)
			prisms_vol += original_elements[k]->getContent();
	}

	INFO("volume of tets: %f", tets_vol);
	INFO("volume of prisms: %f", prisms_vol);

	MeshLib::Mesh split_mesh(refineMesh(*original_mesh));

	std::string out_fname(mesh_out_arg.getValue());
	if (out_fname.empty()) {
		out_fname = BaseLib::dropFileExtension(mesh_out_arg.getValue());
		out_fname += "_refined.vtu";
	}

	FileIO::BoostVtuInterface mesh_io;
	mesh_io.setMesh(&split_mesh);
	mesh_io.writeToFile(out_fname);

	delete original_mesh;
	delete logogCout;
	delete custom_format;
	LOGOG_SHUTDOWN();
}


