/**
 * \file
 * \author Thomas Fischer
 * \date   2010-04-22
 * \brief  Implementation of the Surface class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <list>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "Surface.h"

// GeoLib
#include "AABB.h"
#include "Polygon.h"
#include "SurfaceGrid.h"
#include "AnalyticalGeometry.h"

#include "Triangle.h"
#include "Polyline.h"
#include "EarClippingTriangulation.h"

namespace GeoLib
{
Surface::Surface (const std::vector<Point*> &pnt_vec) :
	GeoObject(), _sfc_pnts(pnt_vec), _bounding_volume(nullptr),
	_surface_grid(nullptr)
{}

Surface::~Surface ()
{
	for (std::size_t k(0); k < _sfc_triangles.size(); k++)
		delete _sfc_triangles[k];
	delete _bounding_volume;
	delete _surface_grid;
}

void Surface::addTriangle(std::size_t pnt_a, std::size_t pnt_b, std::size_t pnt_c)
{
	assert (pnt_a < _sfc_pnts.size() && pnt_b < _sfc_pnts.size() && pnt_c < _sfc_pnts.size());

	// Check if two points of the triangle have identical IDs
	if (pnt_a == pnt_b || pnt_a == pnt_c || pnt_b == pnt_c)
		return;

	_sfc_triangles.push_back(new Triangle(_sfc_pnts, pnt_a, pnt_b, pnt_c));
	if (!_bounding_volume) {
		std::vector<size_t> ids(3);
		ids[0] = pnt_a;
		ids[1] = pnt_b;
		ids[2] = pnt_c;
		_bounding_volume = new AABB<Point>(_sfc_pnts, ids);
		if (_surface_grid == nullptr) {
			_surface_grid = new SurfaceGrid(this);
		}
	} else {
		bool bbx_updated(_bounding_volume->update(*_sfc_pnts[pnt_a]));
		bbx_updated = bbx_updated || _bounding_volume->update(*_sfc_pnts[pnt_b]);
		bbx_updated = bbx_updated || _bounding_volume->update(*_sfc_pnts[pnt_c]);
		if (bbx_updated) {
			delete _surface_grid;
			_surface_grid = new SurfaceGrid(this);
		} else {
			if (! _surface_grid->sortTriangleInGridCells(_sfc_triangles.back())) {
				ERR("Fatal: Could not insert triangle into surface grid. "
					"To keep things consistent, the triangle is removed from "
					"the surface!");
				_sfc_triangles.pop_back();
			}
		}
	}
}

Surface* Surface::createSurface(const Polyline &ply)
{
	if (!ply.isClosed()) {
		WARN("Error in Surface::createSurface() - Polyline is not closed.");
		return nullptr;
	}

	if (ply.getNumberOfPoints() > 2) {
		// create empty surface
		Surface *sfc(new Surface(ply.getPointsVec()));

		Polygon* polygon (new Polygon (ply));
		polygon->computeListOfSimplePolygons ();

		// create surfaces from simple polygons
		const std::list<GeoLib::Polygon*>& list_of_simple_polygons (polygon->getListOfSimplePolygons());
		for (std::list<GeoLib::Polygon*>::const_iterator simple_polygon_it (list_of_simple_polygons.begin());
			simple_polygon_it != list_of_simple_polygons.end(); ++simple_polygon_it) {

			std::list<GeoLib::Triangle> triangles;
			GeoLib::EarClippingTriangulation(*simple_polygon_it, triangles);

			// add Triangles to Surface
			std::list<GeoLib::Triangle>::const_iterator it (triangles.begin());
			while (it != triangles.end()) {
				sfc->addTriangle ((*it)[0], (*it)[1], (*it)[2]);
				++it;
			}
		}
		delete polygon;
		if (sfc->getNTriangles() == 0) {
			WARN("Surface::createSurface(): Triangulation does not contain any triangle.");
			delete sfc;
			return nullptr;
		}
		return sfc;
	} else {
		WARN("Error in Surface::createSurface() - Polyline consists of less than three points and therefore cannot be triangulated.");
		return nullptr;
	}

}

std::size_t Surface::getNTriangles () const
{
	return _sfc_triangles.size();
}

const Triangle* Surface::operator[] (std::size_t i) const
{
	assert (i < _sfc_triangles.size());
	return _sfc_triangles[i];
}

bool Surface::isPntInBoundingVolume(MathLib::Point3d const& pnt) const
{
	return _bounding_volume->containsPoint (pnt);
}

bool Surface::isPntInSfc(MathLib::Point3d const& pnt) const
{
	return _surface_grid->isPointInSurface(
		pnt, std::numeric_limits<double>::epsilon());
}

const Triangle* Surface::findTriangle (MathLib::Point3d const& pnt) const
{
	for (std::size_t k(0); k<_sfc_triangles.size(); k++) {
		if (_sfc_triangles[k]->containsPoint (pnt)) {
			return _sfc_triangles[k];
		}
	}
	return nullptr;
}

} // end namespace
