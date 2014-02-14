/**
 * \date  2014-02-14
 * \brief Program computes the axis aligned boundig box of the given points.
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
#include "Legacy/OGSIOVer4.h"
#include "XmlIO/Boost/BoostXmlGmlInterface.h"

// GeoLib
#include "GEOObjects.h"
#include "AABB.h"

int main(int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Program computes the axis aligned bounding box of given points",
			' ', "0.1");
	TCLAP::MultiArg<std::string> geo_fnames_arg(
		"g",
		"geometry-files",
		"non-empty sequence of file names",
		true,
		"file name strings");
	cmd.add(geo_fnames_arg);
	cmd.parse(argc, argv);

	GeoLib::GEOObjects* geo_objs(new GeoLib::GEOObjects);

	std::vector<std::string> fnames(geo_fnames_arg.getValue());
	std::vector<std::string> unique_names;
	unique_names.resize(fnames.size());
	std::vector<std::string> error_strings;
	for (std::size_t k(0); k<fnames.size(); k++) {
		FileIO::Legacy::readGLIFileV4(
			fnames[k],
			geo_objs,
			unique_names[k],
			error_strings);
	}

	std::string merged_geo_name("AllPoints");
	geo_objs->mergeGeometries(unique_names, merged_geo_name);

	std::vector<GeoLib::Point*> const* all_pnts(geo_objs->getPointVec(merged_geo_name));

	GeoLib::AABB<GeoLib::Point> aabb(all_pnts->begin(), all_pnts->end());
	GeoLib::Point const& min(aabb.getMinPoint());
	GeoLib::Point const& max(aabb.getMaxPoint());
	INFO("AABB: [%f, %f] x [%f, %f] x [%f,%f]",
		min[0], max[0],
		min[1], max[1],
		min[2], max[2]);

	delete geo_objs;

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
