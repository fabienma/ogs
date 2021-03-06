/**
 * \file
 * \author Thomas Fischer
 * \date   2010-04-22
 * \brief  Definition of the AABB class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef AABB_H_
#define AABB_H_

#include <limits>
#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <cassert>
#include <vector>

#include <logog/include/logog.hpp>

#include "MathLib/Point3d.h"
#include "MathLib/MathTools.h"
#include "Point.h"

namespace GeoLib
{
/**
 *
 * \ingroup GeoLib
 *
 * \brief Class AABB is an axis aligned bounding box around a given
 * set of geometric points of (template) type PNT_TYPE.
 * */
template <typename PNT_TYPE = GeoLib::Point>
class AABB
{
public:
	/**
	 * construction of object, initialization the axis aligned bounding box
	 * */
	AABB(std::vector<PNT_TYPE*> const& pnts, std::vector<std::size_t> const& ids)
	{
		assert(! ids.empty());
		init(pnts[ids[0]]);
		for (std::size_t i=1; i<ids.size(); ++i) {
			update(*(pnts[ids[i]]));
		}
	}

	/**
	 * copy constructor.
	 * @param src an axis aligned bounding box
	 */
	AABB(AABB<PNT_TYPE> const& src) :
		_min_pnt(src._min_pnt), _max_pnt(src._max_pnt)
	{}

	/**
	 * Construction of object using input iterators. In contrast to give a vector
	 * this approach is more generic. You can use every (stl) container and
	 * C arrays as input for constructing the object.
	 * @attention{The constructor requires that std::distance(first, last) > 0.}
	 * @param first the input iterator to the initial position in the sequence
	 * @param last the input iterator to the final position in a container, i.e. [first, last).
	 * @attention{The iterator last must be reachable from first.}
	 */
	template <typename InputIterator>
	AABB(InputIterator first, InputIterator last)
	{
		if (std::distance(first,last) <= 0)
		{
			ERR("AABB::AABB(InputIterator first, InputIterator last): first > last");
			std::abort();
		}
		init(*first);
		InputIterator it(first);
		while (it != last) {
			update(*it);
			it++;
		}
	}

	bool update(PNT_TYPE const & pnt)
	{
		bool updated(false);
		for (std::size_t k(0); k<3; k++) {
			if (pnt[k] < _min_pnt[k]) {
				_min_pnt[k] = pnt[k];
				updated = true;
			}
			if (_max_pnt[k] < pnt[k]) {
				_max_pnt[k] = pnt[k];
				updated = true;
			}
		}
		return updated;
	}

	/**
	 * check if point is in the axis aligned bounding box
	 */
	template <typename T>
	bool containsPoint(T const & pnt) const
	{
		if (pnt[0] < _min_pnt[0] || _max_pnt[0] < pnt[0]) return false;
		if (pnt[1] < _min_pnt[1] || _max_pnt[1] < pnt[1]) return false;
		if (pnt[2] < _min_pnt[2] || _max_pnt[2] < pnt[2]) return false;
		return true;
	}

	/**
	 * returns a point that coordinates are minimal for each dimension
	 * for the given point set
	 * @return a point
	 */
	MathLib::Point3d const& getMinPoint() const { return _min_pnt; }

	/**
	 * returns a point that coordinates are maximal for each dimension
	 * within the given point set
	 * @return a point
	 */
	MathLib::Point3d const& getMaxPoint() const { return _max_pnt; }

	/**
	 * Method checks if the given AABB object is contained within the
	 * AABB represented by this object.
	 * @param other_aabb the AABB to test with
	 * @return true if the other AABB is contained in the AABB
	 * represented by this object
	 */
	bool containsAABB(AABB<PNT_TYPE> const& other_aabb) const
	{
		return containsPoint(other_aabb.getMinPoint()) && containsPoint(other_aabb.getMaxPoint());
	}

protected:
	MathLib::Point3d _min_pnt = MathLib::Point3d{std::array<double,3>{{
		std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max()}}};
	MathLib::Point3d _max_pnt = MathLib::Point3d{std::array<double,3>{{
		std::numeric_limits<double>::lowest(),
		std::numeric_limits<double>::lowest(),
		std::numeric_limits<double>::lowest()}}};
private:
	void init(PNT_TYPE const & pnt)
	{
		_min_pnt[0] = _max_pnt[0] = pnt[0];
		_min_pnt[1] = _max_pnt[1] = pnt[1];
		_min_pnt[2] = _max_pnt[2] = pnt[2];
	}
	void init(PNT_TYPE const * pnt)
	{
		init(*pnt);
	}
	void update(PNT_TYPE const * pnt)
	{
		update (*pnt);
	}
};
} // end namespace

#endif /* AABB_H_ */
