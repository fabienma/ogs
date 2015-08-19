/**
 * \file
 * \date   2014-05-08
 * \brief  Definition of class PolygonWithSegmentMarker.
 *
 * \copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef POLYGONWITHSEGMENTMARKER_H_
#define POLYGONWITHSEGMENTMARKER_H_

#include "Polygon.h"

namespace GeoLib {

class PolygonWithSegmentMarker: public GeoLib::Polygon
{
public:
	PolygonWithSegmentMarker(GeoLib::Polyline const& polyline);
	virtual ~PolygonWithSegmentMarker();
	/**
	 * Method marks the segment (default mark is true).
	 * @param seg_num the segment number that should be marked
	 * @param mark_val the value of the flag (true or false)
	 */
	void markSegment(std::size_t seg_num, bool mark_val = true);
	/**
	 * Method returns the value of the mark for the given segment.
	 * @param seg_num segment number
	 * @return either true if the segment is marked or false else
	 */
	bool isSegmentMarked(std::size_t seg_num) const;

	/**
	 * Method calls @see Polyline::addPoint() and initializes the mark of the
	 * corresponding line segment.
	 * @see Polyline::addPoint()
	 */
	virtual void addPoint(std::size_t pnt_id);

	/**
	 * Method calls the @see Polyline::insertPoint() and initializes the inserted line
	 * segment with the same value the previous line segment had.
	 * @see Polyline::insertPoint()
	 */
	virtual void insertPoint(std::size_t pos, std::size_t pnt_id);

private:
	std::vector<bool> _marker;
};

} // end namespace GeoLib

#endif /* POLYGONWITHSEGMENTMARKER_H_ */