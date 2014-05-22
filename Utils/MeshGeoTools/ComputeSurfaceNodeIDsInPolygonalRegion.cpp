/**
 * @file ComputeSurfaceNodeIDsInPolygonalRegion.cpp
 * @date 2014-05-20
 * @brief Computes mesh node ids of mesh nodes within a polygonal region, that resides on the surface.
 *
 * @copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// STL
#include <vector>
#include <algorithm>
#include <string>

// TCLAP
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "LogogSimpleFormatter.h"
#include "StringTools.h"
#include "FileTools.h"

// FileIO
#include "readMeshFromFile.h"
#include "FileIO/XmlIO/Boost/BoostXmlGmlInterface.h"

// GeoLib
#include "PointWithID.h"
#include "Polygon.h"

// MathLib
#include "Vector3.h"

// MeshLib
#include "Mesh.h"
#include "Node.h"
#include "MeshSurfaceExtraction.h"

void writeToFile(std::string const& id_area_fname, std::string const& gli_fname,
	std::vector<std::pair<std::size_t, double>> const& ids_and_areas,
	std::vector<MeshLib::Node*> const& mesh_nodes)
{
	std::ofstream ids_and_area_out(id_area_fname.c_str());
	std::ofstream gli_out(gli_fname.c_str());
	gli_out << "#POINTS\n";
	for (std::size_t k(0); k<ids_and_areas.size(); k++) {
		ids_and_area_out << ids_and_areas[k].first << " " << ids_and_areas[k].second << "\n";
		gli_out << k << " " << GeoLib::Point(mesh_nodes[ids_and_areas[k].first]->getCoords())
			<< " $NAME " << ids_and_areas[k].first << "\n";
	}
	ids_and_area_out.close();
	gli_out << "#STOP";
	gli_out.close();
}

int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Computes ids of mesh nodes that are in a polygonal region and resides on the surface.", ' ', "0.1");
	TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
	                                     "the name of the file containing the input mesh", true,
	                                     "", "file name of input mesh");
	cmd.add(mesh_in);
	TCLAP::ValueArg<std::string> geo_in("g", "gml-file",
	                                     "the name of the gml file containing the polygons", true,
	                                     "", "file name of input geometry");
	cmd.add(geo_in);
	TCLAP::ValueArg<std::string> out_path("o", "output-path",
	                                      "the name of the file the ids will be written to", false,
	                                      "", "file name of id list");
	cmd.add(out_path);
	cmd.parse(argc, argv);

	if (out_path.getValue().empty()) {
		out_path.getValue() = BaseLib::extractPath(geo_in.getValue());
	}

	MeshLib::Mesh* mesh (FileIO::readMeshFromFile(mesh_in.getValue()));
	INFO("Mesh read: %d nodes, %d elements.", mesh->getNNodes(), mesh->getNElements());

	ProjectData project_data;
	FileIO::BoostXmlGmlInterface xml(project_data);
	xml.readFile(geo_in.getValue());
	std::vector<std::string> geo_names;
	project_data.getGEOObjects()->getGeometryNames(geo_names);
	INFO("Geometry read: %d points, %d polylines.",
		project_data.getGEOObjects()->getPointVec(geo_names[0])->size(),
		project_data.getGEOObjects()->getPolylineVec(geo_names[0])->size());

	MathLib::Vector3 const dir(0.0, 0.0, -1.0);
	std::vector<double> areas;
	MeshLib::MeshSurfaceExtraction::getSurfaceAreaForNodes(
		MeshLib::MeshSurfaceExtraction::getMeshSurface(*mesh, dir), areas);

	std::vector<GeoLib::PointWithID*> all_sfc_pnts(
		MeshLib::MeshSurfaceExtraction::getSurfaceNodes(*mesh, dir)
	);

	std::for_each(all_sfc_pnts.begin(), all_sfc_pnts.end(), [](GeoLib::PointWithID* p) { (*p)[2] = 0.0; });

	std::vector<MeshLib::Node*> const& mesh_nodes(mesh->getNodes());
	GeoLib::PolylineVec const* ply_vec(
		project_data.getGEOObjects()->getPolylineVecObj(geo_names[0])
	);
	std::vector<GeoLib::Polyline*> const& plys(*(ply_vec->getVector()));

	for (std::size_t j(0); j<plys.size(); j++) {
		if (! plys[j]->isClosed()) {
			continue;
		}
		std::string polygon_name;
		ply_vec->getNameOfElement(plys[j], polygon_name);
		// create Polygon from Polyline
		GeoLib::Polygon polygon(*(plys[j]));
		// ids of mesh nodes on surface that are within the given polygon
		std::vector<std::pair<std::size_t, double>> ids_and_areas;
		for (std::size_t k(0); k<all_sfc_pnts.size(); k++) {
			GeoLib::PointWithID const& pnt(*(all_sfc_pnts[k]));
			if (polygon.isPntInPolygon(pnt)) {
				ids_and_areas.push_back(std::pair<std::size_t, double>(pnt.getID(), areas[k]));
			}
		}
		INFO("Surface %s contains %d surface nodes.",
			polygon_name.c_str(),
			ids_and_areas.size()
		);
		std::string id_and_area_fname(out_path.getValue() + polygon_name);
		std::string gli_fname(out_path.getValue() + polygon_name);
		if (polygon_name.empty()) {
			id_and_area_fname += BaseLib::number2str(j) + ".txt";
			gli_fname += BaseLib::number2str(j) + ".gli";
		}
		writeToFile(id_and_area_fname, gli_fname, ids_and_areas, mesh_nodes);
	}

	std::for_each(all_sfc_pnts.begin(), all_sfc_pnts.end(), std::default_delete<GeoLib::PointWithID>());
	delete mesh;
	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
