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
#include "FileFinder.h"

// FileIO
#include "FileTools.h"
#include "readMeshFromFile.h"
#include "Legacy/OGSIOVer4.h"
#include "XmlIO/Qt/XmlGmlInterface.h"

#ifdef OGS_BUILD_GUI
// Gui/DataView
#include "Gui/DataView/GEOModels.h"
#else
// GeoLib
#include "GeoLib/GEOObjects.h"
#endif

// OGS
#include "OGS/ProjectData.h"
#include "OGS/FEMCondition.h"

// MeshLib
#include "Mesh.h"


// Utils/MeshEdit/
#include "ExtractMeshNodes.h"

void createSurfaceFromVerticalPolygon(GeoLib::GEOObjects & geo_objs,
		GeoLib::Polygon const* polygon, std::string & sfc_name)
{
	// create copy of polygon points
	std::vector<GeoLib::Point*> *sfc_pnts(new std::vector<GeoLib::Point*>);
	const size_t n_polygon_pnts(polygon->getNumberOfPoints()-1);
	for (size_t k(0); k<n_polygon_pnts; k++) {
		sfc_pnts->push_back(new GeoLib::Point(*(polygon->getPoint(k))));
	}

	// add points to GEOObject object
	geo_objs.addPointVec(sfc_pnts, sfc_name);

	if (sfc_pnts->size() != n_polygon_pnts) {
		geo_objs.removePointVec(sfc_name);
		return;
	}

	// create the surface
	GeoLib::Surface *sfc(new GeoLib::Surface(*sfc_pnts));
	// deploying the special structure of the polygon to create the surface
	size_t id0(0);
	size_t id1(1);
	size_t id2(2);
	sfc->addTriangle(id0, id1, id2);
	id1 = id2;
	id2 = n_polygon_pnts-1;
	sfc->addTriangle(id0, id1, id2);
	for (size_t k(2); k<n_polygon_pnts/2; k++) {
		// triangulate quadrilateral (n_polygon_pnts-k+1, k, k+1, n_polygon_pnts-k)
		id0 = id2;
		id1 = k;
		id2 = k+1;
		sfc->addTriangle(id0, id1, id2);
		id1 = k+1;
		id2 = n_polygon_pnts-k;
		sfc->addTriangle(id0, id1, id2);
	}
	std::vector<GeoLib::Surface*> *sfcs(new std::vector<GeoLib::Surface*>);
	sfcs->push_back(sfc);

	std::map<std::string, std::size_t>* sfc_name_map(new std::map<std::string, std::size_t>);
	sfc_name_map->insert(std::pair<std::string, std::size_t>(sfc_name, 0));

	geo_objs.addSurfaceVec(sfcs, sfc_name, sfc_name_map);
}

int main(int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd(
			"Programme calculates a polygon (that is vertical located) out of a polyline embedded in a mesh. \
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

	TCLAP::ValueArg < std::string > sfc_out_arg("s", "surface-file", "", false, "", "gml file name");
	cmd.add(sfc_out_arg);


	cmd.parse(argc, argv);

	// *** read mesh
	std::string tmp(mesh_arg.getValue());

	std::string file_base_name;
	if (mesh_arg.getValue().find(".msh") != std::string::npos)
		file_base_name = mesh_arg.getValue().substr(0, mesh_arg.getValue().size() - 4);
	MeshLib::Mesh* mesh(FileIO::readMeshFromFile(mesh_arg.getValue()));

	// *** read geometry
	ProjectData project;
#ifdef OGS_BUILD_GUI
	GEOModels* geo(dynamic_cast<GEOModels*>(project.getGEOObjects()));
#else
	GeoLib::GEOObjects* geo(project.getGEOObjects());
#endif
	std::string unique_name;
	std::vector<std::string> error_strings;
	FileIO::Legacy::readGLIFileV4(geo_arg.getValue(), geo, unique_name, error_strings);

	// *** get Polygon
	const std::vector<GeoLib::Polyline*>* plys(geo->getPolylineVec(unique_name));
	INFO("Fetched %d polylines.", plys->size());
	if (!plys) {
		ERR("Could not get vector of polylines.");
		delete mesh;
		return -1;
	}

	MeshLib::ExtractMeshNodes extract_mesh_nodes(mesh);

	// *** generate a polygon from polyline
	std::size_t ply_id_upper(std::min(plys->size(), ply_id_upper_arg.getValue()));
	std::size_t ply_id_lower(ply_id_lower_arg.getValue());

	BaseLib::FileFinder file_finder;
	file_finder.addDirectory(".");
	file_finder.addDirectory("/home/fischeth/workspace/OGS-6/sources/FileIO");


	std::vector<std::string> sfc_names;
	std::string prefix_name (BaseLib::extractBaseNameWithoutExtension(mesh_arg.getValue()));
	for (size_t k(ply_id_lower); k < ply_id_upper; k++) {
		bool closed((*plys)[k]->isClosed());
		if (!closed) {
			std::string sfc_name;
			geo->getPolylineVecObj(unique_name)->getNameOfElement((*plys)[k], sfc_name);
			INFO("Converting polyline %d (%s) to polygon (closed polyline).", k, sfc_name.c_str());

			GeoLib::Polygon* polygon(nullptr);
			extract_mesh_nodes.getPolygonFromPolyline(*((*plys)[k]), geo, unique_name, polygon);

			if (polygon) {
				sfc_name = prefix_name + "-" + sfc_name;
				createSurfaceFromVerticalPolygon(*geo, polygon, sfc_name);

				delete polygon;

				if (geo->getPointVec(sfc_name)) {
					INFO("Adding surface %s.", sfc_name.c_str());
					sfc_names.push_back(sfc_name);
				}
			} else {
				INFO("Surface %s could not be created.", sfc_name.c_str());
			}
		}
	}

	if (! sfc_names.empty()) {
		std::string sfc_project_name(prefix_name + "-Surfaces");
		if (sfc_names.size() > 1)
			geo->mergeGeometries(sfc_names, sfc_project_name);
		else
			sfc_project_name = sfc_names[0];

		std::string schema_name(file_finder.getPath("OpenGeoSysGLI.xsd"));
		FileIO::XmlGmlInterface xml(*geo);
		xml.setNameForExport(sfc_project_name);
		xml.writeToFile(sfc_out_arg.getValue());
	}

	delete mesh;

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
