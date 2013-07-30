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
#include "MeshSurfaceExtraction.h"

// Utils/FileConverter
#include "GocadSGridReader.h"

void markElementsWithFaceSetNodes(MeshLib::Mesh &mesh, std::vector<unsigned> & face_set_prop)
{
	std::size_t const n_elements(mesh.getNElements());
	for (std::size_t k(0); k < mesh.getNNodes(); k++) {
		MeshLib::GocadNode* gocad_node(
				dynamic_cast<MeshLib::GocadNode*>(const_cast<MeshLib::Node*>(mesh.getNode(k))));

		std::size_t const face_set_number(gocad_node->getFaceSetNumber());
		if (face_set_number != std::numeric_limits < std::size_t > ::max()) {
			auto neighbor_elements = gocad_node->getElements();
			for (auto it(neighbor_elements.begin()); it != neighbor_elements.end(); it++) {
				if (*it) {
					std::size_t const element_id((*it)->getValue());
					if (element_id < n_elements)
						face_set_prop[element_id] = face_set_number;
				}
			}
		}
	}
	std::string name("FaceSetElements");
	mesh.addPropertyVec(name, face_set_prop);
}

void generateFaceSetMeshes(MeshLib::Mesh &mesh)
{
	std::size_t const n_elements(mesh.getNElements());
	std::vector<unsigned> face_set_prop(n_elements);
	std::fill(face_set_prop.begin(), face_set_prop.end(), 0);

	markElementsWithFaceSetNodes(mesh, face_set_prop);

	for (std::size_t l(1); l<=2; l++) {
		std::vector<MeshLib::Node*> face_set_nodes;
		std::vector<MeshLib::Element*> face_set_elements;

		std::vector<MeshLib::Node*> remaining_nodes;
		std::vector<MeshLib::Element*> remaining_elements;

		std::vector<MeshLib::Element*> const& elements(mesh.getElements());
		for (std::size_t k(0); k<n_elements; k++) {
			if (face_set_prop[k] != l)
				continue;
			MeshLib::Element const*const elem(elements[k]);
			std::size_t n_faces(elem->getNFaces());
			for (std::size_t j(0); j<n_faces; j++) {
				MeshLib::Element const*const face(elem->getFace(j));
				std::size_t const n_nodes(face->getNNodes());
				std::size_t node_cnt(0); // count nodes belonging to face
				for (std::size_t i(0); i<n_nodes; i++) {
					MeshLib::GocadNode const*const node(
							dynamic_cast<MeshLib::GocadNode*>(
								const_cast<MeshLib::Node*>(face->getNode(i))
							)
					);
					if (node != nullptr) {
						if (node->getFaceSetNumber() != std::numeric_limits<std::size_t>::max())
							node_cnt++;
					}
				}
				if (node_cnt == 4) {
					MeshLib::Element *face_set_elem(face->clone());
					for (std::size_t i(0); i<n_nodes; i++) {
						// deep copy of the face nodes
						face_set_nodes.push_back(new MeshLib::Node(*(face->getNode(i))));
						// reset the node pointer in face_set_elem
						face_set_elem->setNode(i, face_set_nodes[face_set_nodes.size()-1]);
					}
					face_set_elements.push_back(face_set_elem);
				} else {
					if (node_cnt != 2)
						continue;
					MeshLib::Element *remaining_elem(face->clone());
					for (std::size_t i(0); i<n_nodes; i++) {
						// deep copy of the face nodes
						remaining_nodes.push_back(new MeshLib::Node(*(face->getNode(i))));
						// reset the node pointer in remaining_elem
						remaining_elem->setNode(i, remaining_nodes[remaining_nodes.size()-1]);
					}
					remaining_elements.push_back(remaining_elem);
				}
			}
		}

		{
			INFO("Creating face set mesh.");
			MeshLib::Mesh face_set_mesh("GocadSGridFaceSet", face_set_nodes, face_set_elements);
			INFO("Face set mesh created.");

			INFO("Writing face set mesh in vtu format.");
			FileIO::BoostVtuInterface vtu;
			vtu.setMesh(&face_set_mesh);
			// output file name
			std::string mesh_out_fname("FaceSetMesh-" + BaseLib::number2str(l) + ".vtu");
			vtu.writeToFile(mesh_out_fname);
		}

		{
			INFO("Creating remaining mesh.");
			MeshLib::Mesh remaining_mesh("GocadSGridRemaining", remaining_nodes, remaining_elements);
			INFO("remaining mesh created.");

			INFO("Writing remaining mesh in vtu format.");
			FileIO::BoostVtuInterface vtu;
			vtu.setMesh(&remaining_mesh);
			// output file name
			std::string mesh_out_fname("RemainingMesh-" + BaseLib::number2str(l) + ".vtu");
			vtu.writeToFile(mesh_out_fname);
		}
	}
}

void addGocadPropertiesToMesh(FileIO::GocadSGridReader const& reader, MeshLib::Mesh &mesh)
{
	std::vector<std::string> const& prop_names(reader.getPropertyNames());
	for (auto name_it(prop_names.begin()); name_it != prop_names.end(); name_it++) {
		boost::optional<FileIO::GocadSGridReader::GocadProperty const&> prop(reader.getProperty(*name_it));
		if (prop) {
			INFO("Adding property \"%s\".", name_it->c_str());
			mesh.addPropertyVec(*name_it, (*prop)._property_data);
		}
	}
}

MeshLib::Mesh* extractSurfaceMesh(MeshLib::Mesh &mesh)
{
	double* dir (new double[3]);
	dir[0] = 0.0;
	dir[1] = -1.0;
	dir[2] = 0.0;
	return MeshLib::MeshSurfaceExtraction::getMeshSurface(mesh, dir);
}


int main(int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog::Cout *logog(new logog::Cout);
	logog->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Programm reads parts of a Gocad structured mesh", ' ', "0.1");

	// Define a value argument and add it to the command line.
	// A value arg defines a flag and a type of value that it expects
	TCLAP::ValueArg<std::string> sg_file_arg("s", "sg", "structured grid file name", true, "",
			"string");

	// Add the argument sg_file_arg to the CmdLine object. The CmdLine object
	// uses this Arg to parse the command line.
	cmd.add(sg_file_arg);

	cmd.parse(argc, argv);

	// read the Gocad SGrid
	INFO("Start reading Gocad SGrid.");
	FileIO::GocadSGridReader reader(sg_file_arg.getValue());
	INFO("End reading Gocad SGrid.");
	std::vector<MeshLib::Node*> const& nodes(reader.getNodes());
	std::vector<MeshLib::Element*> const& elements(reader.getElements());

	INFO("Creating mesh.");
	MeshLib::Mesh mesh("GocadSGrid", nodes, elements);
	INFO("Mesh created.");

//	INFO("Generating a mesh for every face set.");
//	generateFaceSetMeshes(mesh);
	INFO("Add Gocad properties to mesh.");
	addGocadPropertiesToMesh(reader, mesh);

	{
		MeshLib::Mesh *surface_mesh(extractSurfaceMesh(mesh));
		INFO("Writing surface mesh in vtu format.");
		FileIO::BoostVtuInterface vtu;
		vtu.setMesh(surface_mesh);
		// output file name
		std::string mesh_out_fname(BaseLib::extractPath(sg_file_arg.getValue()) +
					BaseLib::extractBaseNameWithoutExtension(sg_file_arg.getValue()) + "_surface.vtu");
		vtu.writeToFile(mesh_out_fname);
		delete surface_mesh;
	}

	INFO("Writing mesh in vtu format.");
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
