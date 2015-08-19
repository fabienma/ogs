/**
 * \file
 * \date   2014-05-08
 * \brief  Implementation of class PolygonWithSegmentMarker.
 *
 * \copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PolygonWithSegmentMarker.h"

namespace GeoLib {

PolygonWithSegmentMarker::PolygonWithSegmentMarker(GeoLib::Polyline const& polyline)
	: GeoLib::Polygon(polyline, true)
{
	const size_t n_pnts(getNumberOfPoints());
	_marker.resize(n_pnts);
	for (size_t k(0); k<n_pnts; k++) {
		_marker[k] = false;
	}
}

PolygonWithSegmentMarker::~PolygonWithSegmentMarker()
{}


void PolygonWithSegmentMarker::markSegment(size_t seg_num, bool mark_val)
{
	_marker[seg_num] = mark_val;
}
bool PolygonWithSegmentMarker::isSegmentMarked(size_t seg_num) const
{
	return _marker[seg_num];
}

void PolygonWithSegmentMarker::addPoint(size_t pnt_id)
{
	Polyline::addPoint(pnt_id);
	_marker.push_back(false);
}

void PolygonWithSegmentMarker::insertPoint(size_t pos, size_t pnt_id)
{
	Polyline::insertPoint(pos, pnt_id);
	_marker.insert(_marker.begin()+pos, _marker[pos]);
}

} // end GeoLib
