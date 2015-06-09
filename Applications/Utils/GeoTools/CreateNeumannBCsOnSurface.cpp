/*
 * \date 2015-04-27
 * \brief Creates Neumann boundary conditions for all surface mesh nodes.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <cstdlib>
#include <string>
#include <vector>

// TCLAP
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "BaseLib/LogogSimpleFormatter.h"
#include "BaseLib/FileTools.h"

// FileIO
#include "FileIO/Legacy/OGSIOVer4.h"
#include "FileIO/readMeshFromFile.h"
//#include "FileIO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "FileIO/VtkIO/VtuInterface.h"

// GeoLib
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "GeoLib/PointVec.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/convertMeshToGeo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Prism.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Elements/Line.h"

std::string composeName(std::string const& prefix, std::size_t number)
{
	return std::string(prefix+std::to_string(number));
}

int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd(
		"Creates Neumann boundary conditions for all surface nodes",
		' ',
		"0.1");

	TCLAP::ValueArg<std::string> mesh_arg("i", "input-mesh-file",
		"the name of the file containing the mesh", true,
		"", "mesh file name");
	cmd.add(mesh_arg);

	TCLAP::ValueArg<std::string> bc_file_arg("o", "output-bc-file",
		"the name of the output file containing the boundary conditions", true,
		"", "bc file name");
	cmd.add(bc_file_arg);

	TCLAP::ValueArg<std::string> gli_file_arg("g", "geometry-gli-file",
		"the name of the output file containing the boundary conditions", true,
		"", "gli file name");
	cmd.add(gli_file_arg);

	cmd.parse(argc, argv);

	INFO("Reading mesh \"%s\" ... ", mesh_arg.getValue().c_str());
	MeshLib::Mesh * subsfc_mesh(FileIO::readMeshFromFile(mesh_arg.getValue()));
	INFO("done.");

	INFO("Extracting top surface of mesh \"%s\" ... ",
		mesh_arg.getValue().c_str());
	const MathLib::Vector3 dir(0,0,-1);
	double const angle(90);
	std::unique_ptr<MeshLib::Mesh> sfc_mesh(
		MeshLib::MeshSurfaceExtraction::getMeshSurface(
			*subsfc_mesh, dir, angle, true
		)
	);
	INFO("done.");

	// create mapping used for names
	std::vector<std::size_t> sfc_pnt_to_mesh_node_map(sfc_mesh->getNNodes());
	for (std::size_t k(0); k<sfc_mesh->getNNodes(); ++k)
		sfc_pnt_to_mesh_node_map[k] = sfc_mesh->getNodes()[k]->getID();

	sfc_mesh->resetNodeIDs();
	delete subsfc_mesh;

	// transfer the mesh to a geometry the boundary conditions will be set to
	INFO("Convert surface mesh to geometry.");
	GeoLib::GEOObjects geometries;
	MeshLib::convertMeshToGeo(*sfc_mesh, geometries);
//	FileIO::BoostXmlGmlInterface io(geometries);
//	io.setNameForExport(sfc_mesh->getName());
//	io.writeToFile("Surface.gml");

	std::vector<std::string> geometry_names;
	geometries.getGeometryNames(geometry_names);
	std::string const& geometry_name(geometry_names.back());
	INFO("done.");

	// in order to be able to reference the geomtry in the boundary condition
	// set names for all GeoLib::Points (= MeshLib::Nodes)
	INFO("Set names for points in geometry.");
	GeoLib::PointVec & pnt_vec(*(geometries.getPointVecObj(geometry_name)));
	for (std::size_t k(0); k<sfc_mesh->getNNodes(); ++k) {
		std::size_t const id(sfc_pnt_to_mesh_node_map[k]);
		pnt_vec.setNameOfElementByID(k, composeName("N", id));
	}
	INFO("done.");

	// write geometry as a gli-file
	INFO("Write geomtry \"%s\" to file \"%s\".", geometry_name.c_str(),
		gli_file_arg.getValue().c_str());
	FileIO::Legacy::writeGLIFileV4(gli_file_arg.getValue(), geometry_name, geometries);
	INFO("done.");

	// get the area
	INFO("Get the area for the surface nodes.");
	std::vector<double> sfc_area_for_nodes(
		MeshLib::MeshSurfaceExtraction::getSurfaceAreaForNodes(*sfc_mesh));
	INFO("done.");

	// write the Neumann BC to an OGS-5 st file
	INFO("Write the Neumann boundary conditions as OGS-5 source term file "
		"\"%s\".", bc_file_arg.getValue().c_str());
	std::ofstream st_out(bc_file_arg.getValue());
	if (!st_out) {
		ERR("Could not open st file \"%s\" for writing.",
			bc_file_arg.getValue().c_str());
		return EXIT_FAILURE;
	}
	for (std::size_t k(0); k<sfc_mesh->getNNodes(); ++k) {
		std::size_t const id(sfc_pnt_to_mesh_node_map[k]);
		st_out << "#SOURCE_TERM\n";
		st_out << " $PCS_TYPE\n";
		st_out << "  LIQUID_FLOW\n";
		st_out << " $PRIMARY_VARIABLE\n";
		st_out << "  PRESSURE1\n";
		st_out << " $GEO_TYPE\n";
		st_out << "  POINT " << composeName("N", id) << "\n";
		st_out << " $DIS_TYPE CONSTANT\n";
		st_out << "  " << sfc_area_for_nodes[k] << "\n";
	}
	st_out << "#STOP\n";
	st_out.close();

	return EXIT_SUCCESS;
}
