/**
 * \author Thomas Fischer
 * \date   2012-06-05
 * \brief  Programme to convert a open polyline to a vertical polygon
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// STL
#include <iostream>
#include <vector>
#include <string>

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "tclap/CmdLine.h"

// BaseLib
#include "LogogSimpleFormatter.h"

// FileIO
#include "FileTools.h"
#include "readMeshFromFile.h"
#include "Legacy/OGSIOVer4.h"

// GeoLib
#include "GEOObjects.h"

// MeshLib
#include "Mesh.h"

// Utils/MeshEdit/
#include "ExtractMeshNodes.h"

int main(int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd(
			"Programme calculates out of a polyline embedded in a mesh a polygon (that is vertical located). \
			For this reason the polyline is projected to the bottom surface and top surface \
			of the mesh. The resulting two polylines are linked to a closed polyline, i.e., a polygon.",
			' ', "0.1");

	TCLAP::ValueArg<std::string> mesh_arg("m", "mesh",
			"*layered* mesh", true, "", "filename for layered mesh");
	cmd.add(mesh_arg);

	TCLAP::ValueArg<std::string> geo_arg("", "geometry-file", "",
			true, "", "file name");
	cmd.add(geo_arg);

	TCLAP::ValueArg<std::size_t> ply_id_upper_arg("", "upper", "upper boundary of polyline range",
			true, 0, "");
	cmd.add(ply_id_upper_arg);

	TCLAP::ValueArg < std::size_t > ply_id_lower_arg("", "lower", "lower boundary of polyline range", true, 0, "");
	cmd.add(ply_id_lower_arg);

	cmd.parse(argc, argv);

	// *** read mesh
	std::string tmp(mesh_arg.getValue());

	std::string file_base_name;
	if (mesh_arg.getValue().find(".msh") != std::string::npos)
		file_base_name = mesh_arg.getValue().substr(0, mesh_arg.getValue().size() - 4);
	MeshLib::Mesh* mesh(FileIO::readMeshFromFile(mesh_arg.getValue()));

	// *** read geometry
	GeoLib::GEOObjects* geo(new GeoLib::GEOObjects);
	std::string unique_name;
	std::vector<std::string> error_strings;
	FileIO::readGLIFileV4(geo_arg.getValue(), geo, unique_name, error_strings);

	// *** get Polygon
	const std::vector<GeoLib::Polyline*>* plys(geo->getPolylineVec(unique_name));
	std::cout << "fetched polylines: " << std::flush << plys->size() << std::endl;
	if (!plys) {
		std::cout << "could not get vector of polylines" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}

	MeshLib::ExtractMeshNodes extract_mesh_nodes(mesh);

	// *** generate a polygon from polyline
	std::size_t ply_id_upper(std::min(plys->size(), ply_id_upper_arg.getValue()));
	std::size_t ply_id_lower(ply_id_lower_arg.getValue());

	for (size_t k(ply_id_lower); k < ply_id_upper; k++) {
		bool closed((*plys)[k]->isClosed());
		if (!closed) {
			std::cout << "converting polyline " << k << " to polygon (closed polyline) "
							<< std::endl;
			GeoLib::Polygon* polygon(NULL);
			extract_mesh_nodes.getPolygonFromPolyline(*((*plys)[k]), geo, unique_name, polygon);
			std::string *polygon_name(new std::string);
			geo->getPolylineVecObj(unique_name)->getNameOfElementByID(k, *polygon_name);
			(*polygon_name) += "-Polygon";
			geo->getPolylineVecObj(unique_name)->push_back(polygon, polygon_name);
		}
	}

	std::string path(BaseLib::extractPath(mesh_arg.getValue()));
	std::string const fname_of_new_file(path + "New.gli");
	FileIO::writeGLIFileV4(fname_of_new_file, unique_name, *geo);

	delete mesh;
	delete geo;

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
