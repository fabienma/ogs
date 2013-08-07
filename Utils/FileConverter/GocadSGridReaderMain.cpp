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
#include "Elements/Edge.h"
#include "MeshSurfaceExtraction.h"

// Utils/FileConverter
#include "GocadSGridReader.h"

void writeFaceSetNodes(MeshLib::Mesh const& mesh, std::size_t face_set_number, std::string const& path)
{
	std::stringstream ss;
	std::size_t cnt(0); // count face set nodes
	for (std::size_t k(0); k < mesh.getNNodes(); k++) {
		MeshLib::GocadNode* gocad_node(
				dynamic_cast<MeshLib::GocadNode*>(const_cast<MeshLib::Node*>(mesh.getNode(k))));

		bool const face_set_member(gocad_node->isMemberOfFaceSet(face_set_number));
		if (face_set_member) {
			ss << cnt << " " << *gocad_node << " $NAME " << face_set_number << "\n";
			cnt++;
		}
	}
	if (cnt > 0) {
		std::string fname(path + "Surfaces/FaceSetNodes-" + BaseLib::number2str(face_set_number) + ".gli");
		INFO("Writing nodes of face set to file \"%s\".", fname.c_str());
		std::ofstream os(fname.c_str());
		os << "#POINTS\n";
		os << ss.str();
		os << "#STOP";
		os.close();
	}
}

void markElementsWithFaceSetNodes(MeshLib::Mesh const& mesh, std::vector<std::bitset<128>> & face_set_membership)
{
	std::size_t const n_elements(mesh.getNElements());
	for (std::size_t k(0); k < mesh.getNNodes(); k++) {
		MeshLib::GocadNode* gocad_node(
				dynamic_cast<MeshLib::GocadNode*>(const_cast<MeshLib::Node*>(mesh.getNode(k))));

		bool const face_set_member(gocad_node->isMemberOfAnyFaceSet());
		if (face_set_member) {
			auto neighbor_elements = gocad_node->getElements();
			for (auto it(neighbor_elements.begin()); it != neighbor_elements.end(); it++) {
				if (*it) {
					std::size_t const element_id((*it)->getValue());
					if (element_id < n_elements) {
						face_set_membership[element_id]
						              |= gocad_node->getFaceSetMembership();
					}
				}
			}
		}
	}
}

std::size_t getNumberOfNodesInFaceBelongingToFaceSet(MeshLib::Element const* const face,
		std::size_t face_set_number)
{
	std::size_t const n_face_nodes(face->getNNodes());
	std::size_t node_cnt(0); // count nodes belonging to face
	for (std::size_t k(0); k<n_face_nodes; k++) {
		MeshLib::GocadNode const*const node(
				dynamic_cast<MeshLib::GocadNode*>(const_cast<MeshLib::Node*>(face->getNode(k))));
		if (node != nullptr) {
			if (node->isMemberOfFaceSet(face_set_number))
				node_cnt++;
		}
	}

	return node_cnt;
}

void addFaceSetFace(MeshLib::Element const*const face,
		std::vector<MeshLib::Node*> &face_set_nodes,
		std::vector<MeshLib::Element*> &face_set_elements)
{
	MeshLib::Element *face_set_elem(face->clone());
	std::size_t const n_nodes(face_set_elem->getNNodes());
	for (std::size_t i(0); i<n_nodes; i++) {
		// deep copy of the face nodes
		face_set_nodes.push_back(new MeshLib::Node(*(face->getNode(i))));
		// reset the node pointer in face_set_elem
		face_set_elem->setNode(i, face_set_nodes[face_set_nodes.size()-1]);
	}
	face_set_elements.push_back(face_set_elem);
}

bool operator== (MeshLib::Edge const& e0, MeshLib::Edge const& e1)
{
	if (e0.getNode(0) == e1.getNode(0) && e0.getNode(1) == e1.getNode(1))
		return true;
	if (e0.getNode(0) == e1.getNode(1) && e0.getNode(1) == e1.getNode(0))
		return true;

	return false;
}

void generateFaceSetMeshes(MeshLib::Mesh &mesh, std::string const& path)
{
	std::size_t const n_elements(mesh.getNElements());
	std::vector<std::bitset<128>> face_set_membership(n_elements);

	markElementsWithFaceSetNodes(mesh, face_set_membership);

	for (std::size_t l(0); l<128; l++) {
		writeFaceSetNodes(mesh, l, path);
		std::vector<MeshLib::Node*> face_set_nodes;
		std::vector<MeshLib::Element*> face_set_elements;

		std::vector<MeshLib::Element*> const& elements(mesh.getElements());
		for (std::size_t k(0); k<n_elements; k++) {
			if (face_set_membership[k].none())
				continue;
			MeshLib::Element const*const elem(elements[k]);
			std::size_t n_faces(elem->getNFaces());
			for (std::size_t j(0); j<n_faces; j++) {
				MeshLib::Element const*const face(elem->getFace(j));
				std::size_t node_cnt(getNumberOfNodesInFaceBelongingToFaceSet(face, l));

				if (node_cnt == 4) {
					addFaceSetFace(face, face_set_nodes, face_set_elements);
					continue;
				}

				// check the neighbor to fix gaps
				if (node_cnt == 2) {
					std::size_t const n_neighbor_elems(elem->getNNeighbors());
					for (std::size_t i(0); i<n_neighbor_elems; i++) {
						MeshLib::Element const*const neighbor_elem(elem->getNeighbor(i));
						if (neighbor_elem == nullptr)
							continue;
						std::size_t n_neighbor_faces(neighbor_elem->getNFaces());
						for (std::size_t jj(0); jj < n_neighbor_faces; jj++) {
							MeshLib::Element const* const neighbor_face(neighbor_elem->getFace(jj));
							std::size_t const neighbor_node_cnt(getNumberOfNodesInFaceBelongingToFaceSet(neighbor_face, l));
							if (neighbor_node_cnt >= 2) {
								addFaceSetFace(face, face_set_nodes, face_set_elements);
								addFaceSetFace(neighbor_face, face_set_nodes, face_set_elements);
							}
						}
					}
				}
			}
		}

		{
			if (face_set_nodes.size() == 0)
				continue;
			INFO("Creating face set mesh.");
			MeshLib::Mesh face_set_mesh("GocadSGridFaceSet", face_set_nodes, face_set_elements);
			INFO("Face set mesh created. #nodes: %d, #elements: %d", face_set_mesh.getNNodes(),
					face_set_mesh.getNElements());

			FileIO::BoostVtuInterface vtu;
			vtu.setMesh(&face_set_mesh);
			// output file name
			std::string mesh_out_fname(path+"Surfaces/FaceSetMesh-" + BaseLib::number2str(l) + ".vtu");
			INFO("Writing face set mesh \"%s\" in vtu format.", mesh_out_fname.c_str());
			vtu.writeToFile(mesh_out_fname);
		}
	}
}

void cleanUpNoDataValues(MeshLib::Mesh &mesh, double no_data_value,
		std::vector<double> const& original,
		std::vector<double> & reworked)
{
	std::vector<MeshLib::Element*> const& elements(mesh.getElements());
	for (std::size_t k(0); k<original.size(); ++k) {
		if (std::abs(original[k] - no_data_value) > std::numeric_limits<double>::epsilon()) {
			reworked[k] = original[k];
		} else {
			std::size_t const n_neighbors(elements[k]->getNNeighbors());
			double prop_value(0.0);
			std::size_t cnt(0); // count neighbors with valid property values
			// simple average
			for (std::size_t j(0); j<n_neighbors; ++j) {
				double neighbor_val(original[elements[k]->getNeighbor(j)->getValue()]);
				if (std::abs(neighbor_val - no_data_value) > std::numeric_limits<double>::epsilon()) {
					prop_value += neighbor_val;
					cnt++;
				}
			}
			prop_value /= cnt;
			reworked[k] = prop_value;
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
			std::vector<double> reworked_properties;
			reworked_properties.resize((*prop)._property_data.size());
			cleanUpNoDataValues(mesh, (*prop)._property_no_data_value, (*prop)._property_data, reworked_properties);
			mesh.addPropertyVec(*name_it, reworked_properties);
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

void writeMeshPropertiesToFile(std::string const& fname, std::string const& prop_name,
		MeshLib::Mesh const& mesh)
{
	boost::optional<std::vector<double> const&> prop_vec(mesh.getDoublePropertyVec(prop_name));
	if (! prop_vec) {
		ERR("Could not read property \"%s\" from mesh.", prop_name.c_str());
		return;
	}

	std::ofstream os(fname.c_str());
	if (!os) {
		ERR("Could not open file \"%s\".", fname.c_str());
		return;
	} else {
		INFO("Open file \"%s\" for writing.", fname.c_str());
	}

	os << "# OpenGeoSys material property file - " << prop_name << " \n";
	if (prop_name.compare("Porosity") == 0) {
		for (std::size_t k(0); k<(*prop_vec).size(); k++) {
			os << k << " " << (*prop_vec)[k] / 100.0 << "\n";
		}
	}
	if (prop_name.compare("Permeability") == 0) {
		const double scale(9.86923 * 1e-16 * 9.81 * 1e3 * 1e3);
		for (std::size_t k(0); k<(*prop_vec).size(); k++) {
			os << k << " " << (*prop_vec)[k] * scale << "\n";
		}
	}

	os << "#STOP\n";
	os.close();
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

	INFO("Add Gocad properties to mesh.");
	addGocadPropertiesToMesh(reader, mesh);

	INFO("Generating a mesh for every face set.");
	generateFaceSetMeshes(mesh, BaseLib::extractPath(sg_file_arg.getValue()));

//	{
//		MeshLib::Mesh *surface_mesh(extractSurfaceMesh(mesh));
//		INFO("Writing surface mesh in vtu format.");
//		FileIO::BoostVtuInterface vtu;
//		vtu.setMesh(surface_mesh);
//		// output file name
//		std::string mesh_out_fname(BaseLib::dropFileExtension(sg_file_arg.getValue()) + "_surface.vtu");
//		vtu.writeToFile(mesh_out_fname);
//		delete surface_mesh;
//	}

	{
		std::string prop_name("Porosity");
		std::string prop_fname(BaseLib::dropFileExtension(sg_file_arg.getValue()) +
						"-" + prop_name + ".txt");
		writeMeshPropertiesToFile(prop_fname, prop_name, mesh);
	}

	{
		std::string prop_name("Permeability");
		std::string prop_fname(BaseLib::dropFileExtension(sg_file_arg.getValue()) +
						"-" + prop_name + ".txt");
		writeMeshPropertiesToFile(prop_fname, prop_name, mesh);
	}

	INFO("Writing mesh in vtu format.");
	FileIO::BoostVtuInterface vtu;
	vtu.setMesh(&mesh);
	// output file name
	std::string mesh_out_fname(BaseLib::dropFileExtension(sg_file_arg.getValue()) + ".vtu");

	vtu.writeToFile(mesh_out_fname);

	delete logog;
	delete custom_format;
	LOGOG_SHUTDOWN();

	return 0;
}
