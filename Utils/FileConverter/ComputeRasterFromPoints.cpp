/**
 * \date  2014-02-14
 * \brief Program computes a raster out of a given point set.
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

// GeoLib
#include "GEOObjects.h"
#include "Raster.h"

int main(int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd(
		"Program computes a raster out of a given point set",
		' ',
		"0.1");
	TCLAP::ValueArg<std::string> geo_fname_arg(
		"g",
		"geometry-file",
		"non-empty file name",
		true,
		"",
		"file name string");
	cmd.add(geo_fname_arg);
	TCLAP::ValueArg<unsigned> n_rows_arg(
		"r",
		"nRows",
		"number of rows",
		true,
		10u,
		"integer number used as number of rows");
	cmd.add(n_rows_arg);
	TCLAP::ValueArg<unsigned> n_cols_arg(
		"c",
		"nCols",
		"number of columns",
		true,
		10u,
		"integer number used as number of columns");
	cmd.add(n_cols_arg);

	TCLAP::ValueArg<double> x_min_arg(
		"",
		"xmin",
		"minimal value of x ",
		false,
		0.0,
		"floating point number for the minimal value for x");
	cmd.add(x_min_arg);
	TCLAP::ValueArg<double> y_min_arg(
		"",
		"ymin",
		"minimal value of y",
		false,
		0.0,
		"floating point number for the minimal value for y");
	cmd.add(y_min_arg);
	TCLAP::ValueArg<double> edge_size_of_cell_arg(
		"e",
		"edge-size",
		"size of the edge of a cell",
		false,
		1.0,
		"floating point number for the size of an edge of a cell");
	cmd.add(edge_size_of_cell_arg);
	TCLAP::ValueArg<std::size_t> dilate_arg(
		"d",
		"dilate-pixels",
		"number of neighbor pixels that have no-data values and that have a neighbor neighbor pixel with a valid value",
		false,
		0,
		"positiv integer number (of neighbor pixels to dilate)");
	cmd.add(dilate_arg);


	cmd.parse(argc, argv);

	GeoLib::GEOObjects* geo_objs(new GeoLib::GEOObjects);

	std::string fname(geo_fname_arg.getValue());
	std::string unique_name;
	std::vector<std::string> error_strings;
	FileIO::Legacy::readGLIFileV4(
		fname,
		geo_objs,
		unique_name,
		error_strings);

	// init raster data
	const double no_data_val(-9999.0);
	std::vector<double> data;
	const std::size_t n_rows(n_rows_arg.getValue());
	const std::size_t n_cols(n_cols_arg.getValue());
	data.resize(n_rows * n_cols);
	std::fill(data.begin(), data.end(), no_data_val);

	std::vector<GeoLib::Point*> const* pnts(geo_objs->getPointVec(unique_name));
	const double x_min(x_min_arg.getValue());
	const double y_min(y_min_arg.getValue());
	const double cell_size(edge_size_of_cell_arg.getValue());
	for (std::size_t k(0); k<pnts->size(); k++) {
		GeoLib::Point const& pnt(*(*pnts)[k]);
		// compute pos in grid
		const double x((pnt[0] - x_min)/cell_size - 0.5);
		const double y((pnt[1] - y_min)/cell_size - 0.5);
		if (x<0 || y<0 || x>n_cols || y>n_rows)
			continue;
		const std::size_t xpix(static_cast<std::size_t>(x));
		const std::size_t ypix(static_cast<std::size_t>(y));
		data[xpix + ypix*n_cols] = pnt[2];
	}

	const double eps (std::numeric_limits<double>::epsilon());
	for (std::size_t k(0); k<dilate_arg.getValue(); k++) {
		// dilate one layer
		for (std::size_t i(1); i<n_cols; i++) {
			for (std::size_t j(1); j<n_rows; j++) {
				if (data[j*n_cols+i] != no_data_val) {
					if (fabs(data[j*n_cols+i-1] - no_data_val) < eps)
						data[j*n_cols+i-1] = data[j*n_cols+i];
					if (fabs(data[(j-1)*n_cols+i] - no_data_val) < eps)
						data[(j-1)*n_cols+i] = data[j*n_cols+i];
				}
			}
		}
		for (std::size_t i(0); i<n_cols-1; i++) {
			for (std::size_t j(0); j<n_rows-1; j++) {
				if (data[(n_rows-j-1)*n_cols+(n_cols-i-1)] != no_data_val) {
					if (fabs(data[(n_rows-j)*n_cols+(n_cols-i-1)] - no_data_val) < eps)
						data[(n_rows-j)*n_cols+(n_cols-i-1)] = data[(n_rows-j-1)*n_cols+n_cols-i-1];
				}
			}
		}
	}

	GeoLib::Raster raster(n_cols, n_rows, x_min, y_min, cell_size,
		data.begin(), data.end());

	std::string out_fname(geo_fname_arg.getValue()+".asc");
	std::ofstream os(out_fname.c_str());
	raster.writeRasterAsASC(os);

	delete geo_objs;

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
