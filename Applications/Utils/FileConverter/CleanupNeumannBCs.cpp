/**
 * \file
 * \brief  Util to cleanup the Neumann BC (OGS-5: st file)
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// STL
#include <string>
#include <vector>

// boost
#include <boost/algorithm/string.hpp> // for trim

// ThirdParty
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "BaseLib/FileTools.h"
#include "BaseLib/LogogSimpleFormatter.h"

// FileIO
#include "FileIO/Legacy/OGSIOVer4.h"

// GeoLib
#include "GeoLib/GEOObjects.h"

struct SimpleOGS5SourceTerm
{
	std::string _pcs_type;
	std::string _primary_variable;
	std::string _pnt_name;
	double _value;
};

void write(std::ostream &os, SimpleOGS5SourceTerm const*const st)
{
	os << "#SOURCE_TERM\n";
	os << " $PCS_TYPE\n";
	os << "  LIQUID_FLOW\n";
	os << " $PRIMARY_VARIABLE\n";
	os << "  PRESSURE1\n";
	os << " $GEO_TYPE\n";
	os << "  POINT " << st->_pnt_name << "\n";
	os << " $DIS_TYPE\n";
	os << "  CONSTANT " << st->_value << "\n";
	os << " $TIM_TYPE\n";
	os << "  CURVE 1\n";
}

SimpleOGS5SourceTerm* read(std::istream &in)
{
	std::string line;
	std::getline(in, line);
	boost::trim(line);
	if (line.compare("#STOP") == 0) {
		return nullptr;
	}
	if (line.compare("#SOURCE_TERM") != 0) {
		ERR("Reading Neumann BC \"#SOURCE_TERM\" keyword.");
		return nullptr;
	}

	std::getline(in, line);
	boost::trim(line);
	if (line.compare("$PCS_TYPE") != 0) {
		ERR("Reading Neumann BC: \"$PCS_TYPE\" keyword.");
		return nullptr;
	}

	SimpleOGS5SourceTerm *st(new SimpleOGS5SourceTerm);

	std::getline(in, line);
	boost::trim(line);
	st->_pcs_type = line;

	std::getline(in, line);
	boost::trim(line);
	if (line.compare("$PRIMARY_VARIABLE") != 0) {
		ERR("Reading Neumann BC: \"$PRIMARY_VARIABLE\" keyword.");
		delete st;
		return nullptr;
	}

	std::getline(in, line);
	boost::trim(line);
	st->_pcs_type = line;

	std::getline(in, line);
	boost::trim(line);
	if (line.compare("$GEO_TYPE") != 0) {
		ERR("Reading Neumann BC: \"$GEO_TYPE\" keyword.");
		delete st;
		return nullptr;
	}

	std::getline(in, line);
	std::size_t pos(line.rfind(" ")+1);
	st->_pnt_name = line.substr(pos);

	std::getline(in, line);
	boost::trim(line);
	if (line.compare("$DIS_TYPE CONSTANT") != 0) {
		ERR("Reading Neumann BC: \"$DIS_TYPE CONSTANT\" keyword.");
		delete st;
		return nullptr;
	}

	std::getline(in, line);
	std::stringstream ss(line);
	if (!(ss >> st->_value)) {
		ERR("Reading Neumann BC: value.");
		delete st;
		return nullptr;
	}

	return st;
}

int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *const custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Cleanup the Neumann BCs not found in the geometry file", ' ', "0.1");

	TCLAP::ValueArg<std::string> neumann_out_arg(
		"o",
		"out",
		"filename for output",
		true,
		"",
		"filename as string");
	cmd.add(neumann_out_arg);

	TCLAP::ValueArg<std::string> neumann_in_arg(
		"i",
		"in",
		"file containing Neumann BCs (OGS-5: st-file)",
		true,
		"",
		"filename as string");
	cmd.add(neumann_in_arg);

	TCLAP::ValueArg<std::string> geom_arg(
		"g",
		"geometry",
		"file containing the geometry the Neumann BCs are linked with",
		true,
		"",
		"filename as string");
	cmd.add(geom_arg);

	cmd.parse(argc, argv);

	// *** Neumann BCs
	std::vector<SimpleOGS5SourceTerm*> st_vec;
	INFO("Reading %s.", neumann_in_arg.getValue().c_str());
	std::ifstream in(neumann_in_arg.getValue());
	while (in) {
		SimpleOGS5SourceTerm* st(read(in));
		if (st)
			st_vec.push_back(st);
		else
			ERR("Reading data record %d.", st_vec.size());
	}
	in.close();
	INFO("\tDone, Read %d data records.", st_vec.size());

	// *** geometry
	GeoLib::GEOObjects geom_obj;
	std::vector<std::string> errors;
	std::string unique_name;
	INFO("Reading geometry from \"%s\".", geom_arg.getValue().c_str());
	if (FileIO::Legacy::readGLIFileV4(geom_arg.getValue(), &geom_obj,
		unique_name, errors)) {
		WARN("There were some issues while reading the geometry file.");
	}
	INFO("\tDone.");

	// cleanup
	std::ofstream out(neumann_out_arg.getValue());
	if (!out) {
		std::for_each(st_vec.begin(), st_vec.end(),
			std::default_delete<SimpleOGS5SourceTerm>());
		ERR("Could not open output file \"%s\".",
			neumann_out_arg.getValue().c_str());
		return EXIT_FAILURE;
	}
	GeoLib::PointVec const& pnt_vec(*(geom_obj.getPointVecObj(unique_name)));
	for (auto st : st_vec) {
		if (pnt_vec.getElementByName(st->_pnt_name))
			write(out, st);
		else {
			INFO("geometrical object \"%s\" for neumann condition not found.",
				st->_pnt_name.c_str());
		}
	}
	out << "#STOP\n";
	out.close();

	std::for_each(st_vec.begin(), st_vec.end(),
		std::default_delete<SimpleOGS5SourceTerm>());

	return EXIT_SUCCESS;
}

