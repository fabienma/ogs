/**
 * \brief  Implements algorithms for extracting information from a mesh.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EXTRACTMESHNODES_H_
#define EXTRACTMESHNODES_H_

// STL
#include <iostream>

#include "MeshLib/Mesh.h"

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Polygon.h"

namespace MeshLib
{
/**
 * This class implements an algorithm to extract mesh node ids from a given *extruded* mesh.
 */
class ExtractMeshNodes
{
public:
	/**
	 * constructor - takes a *extruded* mesh
	 * @param msh an instance of class CFEMesh
	 */
	ExtractMeshNodes(const MeshLib::Mesh* mesh);

	/**
	 * Method computes the polygon to a given polyline that is consisting of the projection
	 * of this polyline to the bottom and the top surface of the mesh and the links between
	 * these two polylines.
	 * @param polyline the ("defining") polyline
	 * @param geo_obj geometric objects manager
	 * @param name the name of the group of geometric objects
	 * @param polygon pointer to the resulting polygon
	 *      warning: the pointer to an already existing polygon will be destroyed
	 */
	void getPolygonFromPolyline(const GeoLib::Polyline& polyline,
		GeoLib::GEOObjects & geo_obj, std::string const& name,
		GeoLib::Polygon* &polygon) const;

private:
	/**
	 * computes the mesh nodes along a polyline belonging to the bottom surface
	 * @param ply computation along the polyline ply
	 * @param bottom_points the bottom mesh nodes as points
	 */
	void getBottomMeshNodesAlongPolylineAsPoints(const GeoLib::Polyline& ply,
		std::vector<GeoLib::Point*>& bottom_points) const;

	/**
	 * computes the mesh nodes along a polyline belonging to the top surface
	 * @param polyline computation along the polyline ply
	 * @param top_points the top mesh nodes as points
	 */
	void getTopMeshNodesAlongPolylineAsPoints(
	        const GeoLib::Polyline& polyline,
	        std::vector<GeoLib::Point*>& top_points) const;
	/**
	 * This method searchs all mesh nodes with the same x and y coordinates
	 * like the polyline points. It returns the mesh nodes as points and the ids
	 * within objects of type Point. The mesh nodes / points are
	 * sorted lexicographically.
	 * @param polyline the "defining" polyline
	 * @param nodes_as_points vector of GeoLib::Point objects
	 */
	void getOrthogonalProjectedMeshNodesAlongPolyline (
	        GeoLib::Polyline const& polyline,
	        std::vector<GeoLib::Point>& nodes_as_points) const;

	const MeshLib::Mesh* _mesh;
};

}

#endif /* EXTRACTMESHNODES_H_ */
