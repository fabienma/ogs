/**
 * @date 2014-10-15
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "gtest/gtest.h"

#include "FileIO/TINInterface.h"
#include "GeoLib/Surface.h"

TEST(GeoLib, IsPointOnSurface)
{
	std::vector<GeoLib::Point*> pnt_vec;
	GeoLib::Surface* sfc(
		FileIO::TINInterface::readTIN(
			"/home/fischeth/Documents/visdata/tom/Pipe3D/SURF_AUSSEN.tin",
			pnt_vec
		)
	);

	GeoLib::Point search_pnt(0.920699085827288, 0.651032443723371, 0);
	EXPECT_FALSE(sfc->isPntInSfc(search_pnt));
	search_pnt[0] = 0.0;
	search_pnt[1] = 2.0;
	search_pnt[2] = 0.5;
	EXPECT_TRUE(sfc->isPntInSfc(search_pnt));
}
