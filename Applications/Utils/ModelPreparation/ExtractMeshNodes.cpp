/*
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <algorithm>

#include "ExtractMeshNodes.h"

#include "BaseLib/quicksort.h"

#include "GeoLib/Point.h"

#include "MathLib/MathTools.h"

#include "MeshLib/Node.h"

namespace MeshLib
{
ExtractMeshNodes::ExtractMeshNodes(const MeshLib::Mesh* mesh) :
	_mesh (mesh)
{
}

void ExtractMeshNodes::getOrthogonalProjectedMeshNodesAlongPolyline(
        GeoLib::Polyline const& polyline,
        std::vector<GeoLib::Point>& nodes_as_points) const
{
	// get all nodes of mesh
	std::vector<MeshLib::Node*> const& nodes(_mesh->getNodes());
	std::size_t number_of_mesh_nodes (nodes.size());

	std::size_t number_of_ply_pnts(polyline.getNumberOfPoints());
	if (polyline.isClosed())
		number_of_ply_pnts--;

	const double eps(50); // _mesh->getSearchLength()); ToDo

	for (std::size_t k(0); k < number_of_ply_pnts; k++) {
		GeoLib::Point proj_ply_pnt((*(polyline.getPoint(k)))[0], (*(polyline.getPoint(k)))[1], 0.0);
		for (std::size_t j(0); j < number_of_mesh_nodes; j++) {
			MathLib::Point3d mesh_pnt(*nodes[j]);
			mesh_pnt[2] = 0.0;
			if (MathLib::sqrDist(mesh_pnt, proj_ply_pnt) < eps)
				nodes_as_points.push_back(GeoLib::Point(*nodes[j], j));
		}
	}

	std::vector<std::size_t> perm;
	for (std::size_t k(0); k < nodes_as_points.size(); k++)
		perm.push_back(k);
	BaseLib::quicksort<GeoLib::Point> (nodes_as_points, 0, nodes_as_points.size(), perm);
}

void ExtractMeshNodes::getTopMeshNodesAlongPolylineAsPoints(
        const GeoLib::Polyline& polyline,
        std::vector<GeoLib::Point*>& top_points) const
{
	std::vector<GeoLib::Point> nodes_as_points;
	this->getOrthogonalProjectedMeshNodesAlongPolyline(polyline, nodes_as_points);

	if (nodes_as_points.empty())
		return;

	double eps (std::numeric_limits<double>::epsilon());
	// collect data (lowest points with same x and y coordinates)
	std::size_t upper_bound (nodes_as_points.size() - 1);
	for (std::size_t k(0); k < upper_bound; k++)
	{
		const GeoLib::Point& p0 (nodes_as_points[k]);
		const GeoLib::Point& p1 (nodes_as_points[k + 1]);
		if (fabs (p0[0] - p1[0]) > eps || fabs (p0[1] - p1[1]) > eps)
			top_points.push_back(new GeoLib::Point(nodes_as_points[k]));
	}
	top_points.push_back(new GeoLib::Point(nodes_as_points[upper_bound]));
}

void ExtractMeshNodes::getBottomMeshNodesAlongPolylineAsPoints(
        const GeoLib::Polyline& polyline,
        std::vector<GeoLib::Point*>& bottom_points) const
{
	std::vector<GeoLib::Point> nodes_as_points;
	this->getOrthogonalProjectedMeshNodesAlongPolyline(polyline, nodes_as_points);

	if (nodes_as_points.empty())
		return;

	double eps(std::numeric_limits<double>::epsilon());
	// collect data (lowest points with same x and y coordinates)
	bottom_points.push_back(new GeoLib::Point(nodes_as_points[0]));
	std::size_t upper_bound (nodes_as_points.size() - 1);
	for (std::size_t k(0); k < upper_bound; k++)
	{
		const GeoLib::Point& p0 (nodes_as_points[k]);
		const GeoLib::Point& p1 (nodes_as_points[k + 1]);
		if (fabs (p0[0] - p1[0]) > eps || fabs (p0[1] - p1[1]) > eps)
			bottom_points.push_back(new GeoLib::Point(nodes_as_points[k+1]));
	}
}

void ExtractMeshNodes::getPolygonFromPolyline (const GeoLib::Polyline& polyline,
                                               GeoLib::GEOObjects & geo_obj,
                                               std::string const& name,
                                               GeoLib::Polygon* & polygon) const
{
	std::vector<GeoLib::Point*> top_polygon_pnts;
	this->getTopMeshNodesAlongPolylineAsPoints(polyline, top_polygon_pnts);

	std::vector<GeoLib::Point*> bottom_polygon_pnts;
	this->getBottomMeshNodesAlongPolylineAsPoints(polyline, bottom_polygon_pnts);

	if (top_polygon_pnts.size() != bottom_polygon_pnts.size() || top_polygon_pnts.size() <= 1) {
		return;
	}

	// append new points to the end of the points vector
	std::vector<std::size_t> top_ids;
	GeoLib::PointVec & pnt_vec(*(geo_obj.getPointVecObj(name)));
	for (auto p : top_polygon_pnts) {
		top_ids.push_back(pnt_vec.push_back(p));
	}
	std::vector<std::size_t> bottom_ids;
	for (auto p : bottom_polygon_pnts) {
		bottom_ids.push_back(pnt_vec.push_back(p));
	}

	// create (an empty) polygon
	GeoLib::Polyline ply(*(geo_obj.getPointVec(name)));

	std::vector<GeoLib::Point*> const* orig_pnts(pnt_vec.getVector());

	// *** add ids of points to polygon
	// for top polyline sort points along polyline
	const double eps(50); // _mesh->getSearchLength()); ToDo
	std::size_t s (top_ids.size());
	for (std::size_t j(0); j < polyline.getNumberOfPoints(); j++)
		for (std::size_t k(0); k < s; k++)
		{
			GeoLib::Point test_pnt(*(*orig_pnts)[top_ids[k]]);
			test_pnt[2] = (*polyline.getPoint(j))[2];
			if (MathLib::sqrDist(*polyline.getPoint(j),test_pnt) < eps) {
				ply.addPoint(top_ids[k]);
				k = s;
			}
		}

	// for bottom polyline sort points along polyline in reverse order
	s = bottom_ids.size();
	GeoLib::Point test_pnt(0.0, 0.0, 0.0);
	for (int j(polyline.getNumberOfPoints() - 1); j > -1; j--) {
		for (std::size_t k(0); k < s; k++) {
			test_pnt[0] = (*(*orig_pnts)[bottom_ids[k]])[0];
			test_pnt[1] = (*(*orig_pnts)[bottom_ids[k]])[1];
			test_pnt[2] = (*polyline.getPoint(j))[2];
			if (MathLib::sqrDist(*polyline.getPoint(j), test_pnt) < eps) {
				ply.addPoint(bottom_ids[k]);
				k = s;
			}
		}
	}

	// close polygon
	polygon = new GeoLib::Polygon(ply);
	polygon->addPoint(polygon->getPointID(0));
	polygon->initialise();
}

} // end namespace MeshLib
