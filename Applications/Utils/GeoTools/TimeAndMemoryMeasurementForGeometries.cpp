/**
 * \file   TimeAndMemoryMeasurementForGeometries.cpp
 * \date   2015-05-26
 * \brief  Generates a geometry for time and memory measurements.
 *
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 */

#include <string>
#include <random>
#include <vector>

#include "tclap/CmdLine.h"
#include "logog/include/logog.hpp"

#include "BaseLib/LogogSimpleFormatter.h"
#include "BaseLib/BuildInfo.h"
#include "BaseLib/MemWatch.h"

#include "GeoLib/GEOObjects.h"

int main(int argc, char *argv[])
{
	LOGOG_INITIALIZE();
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog::Cout *logogCout(new logog::Cout);
	logogCout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Generates a geometry and measures time and memory for "
		"data structure initialization", ' ', BaseLib::BuildInfo::git_describe);
	TCLAP::ValueArg<std::size_t> size_arg("p", "number-of-points",
		"number of points", false, 10000, "positive integer number");
	TCLAP::ValueArg<std::size_t> name_arg("n", "number-of-names", "number of "
		"named points", false, 1000, "positive integer number");
	TCLAP::SwitchArg quiet_arg("q", "quiet", "if quiet arg is given only the "
		"memory consumption is outputed", false);
	cmd.add(size_arg);
	cmd.add(name_arg);
	cmd.add(quiet_arg);
	cmd.parse(argc, argv);

	BaseLib::MemWatch mem_watch;
	unsigned long virt_mem_begin(mem_watch.getVirtMemUsage());
	unsigned long res_mem_begin(mem_watch.getResMemUsage());
	unsigned long shr_mem_begin(mem_watch.getShrMemUsage());

	// generate points with random coordinates
	std::size_t const n_pnts(size_arg.getValue());
	std::vector<GeoLib::Point*> *pnts(new std::vector<GeoLib::Point*>(n_pnts));
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, 10000);
	for (std::size_t k(0); k<n_pnts; ++k) {
		(*pnts)[k] = new GeoLib::Point(dis(gen), dis(gen), dis(gen), k);
	}

	// generate names
	std::map<std::string, std::size_t> *pnt_names(
		new std::map<std::string,std::size_t>);
	std::size_t const n_names(name_arg.getValue());
	for (std::size_t k(0); k<n_names; ++k) {
		(*pnt_names)["P"+std::to_string(k)] = (*pnts)[k]->getID();
	}

	{
		if (! quiet_arg.getValue())
			INFO("Number of points: %d, number of names: %d.", n_pnts, n_names);
		GeoLib::GEOObjects geo_objects;
		std::string geo_name("TestPoints");
		geo_objects.addPointVec(pnts, geo_name, pnt_names);

		unsigned long virt_mem_end(mem_watch.getVirtMemUsage());
		unsigned long res_mem_end(mem_watch.getResMemUsage());
		unsigned long shr_mem_end(mem_watch.getShrMemUsage());

		unsigned long const virt_mem(virt_mem_end-virt_mem_begin);
		unsigned long const res_mem(res_mem_end-res_mem_begin);
		unsigned long const shr_mem(shr_mem_end-shr_mem_begin);
		if (! quiet_arg.getValue()) {
			INFO("Memory consumption after GEOObjects.addPointVec() is: %i kB.",
				(virt_mem+res_mem+shr_mem)/1024);
			INFO("\tvirt mem: %i kB, res mem: %i kB, shared mem: %i kB",
				virt_mem/1024, res_mem/1024, shr_mem/1024);
		} else {
			INFO("%i kB", (virt_mem+res_mem+shr_mem)/1024);
		}
	}
	delete logogCout;
	delete custom_format;
	LOGOG_SHUTDOWN();
	return 0;
}
