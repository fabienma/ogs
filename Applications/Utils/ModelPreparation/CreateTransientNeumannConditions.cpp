/**
 * \date   2015-06-25
 * \brief  Creating transient Neumann conditions
 *
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 */

#include <string>
#include <memory>

#include "boost/date_time/gregorian/gregorian.hpp"

#include "logog/include/logog.hpp"
#include "tclap/CmdLine.h"
#include "tclap/MultiArg.h"

#include "BaseLib/BuildInfo.h"
#include "BaseLib/LogogSimpleFormatter.h"
#include "BaseLib/FileTools.h"

#include "FileIO/readMeshFromFile.h"
#include "FileIO/AsciiRasterInterface.h"
#include "FileIO/XmlIO/Boost/BoostXmlGmlInterface.h"

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Raster.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSurfaceExtraction.h"

void writeGeometry(MeshLib::Mesh const& mesh,
	std::vector<std::pair<std::size_t,double>> const& id_val_pairs,
	std::string const& out_fname)
{
	std::vector<GeoLib::Point*> *pnts(new std::vector<GeoLib::Point*>);
	std::map<std::string, std::size_t> *name_id_map(
		new std::map<std::string, std::size_t>);
	for (auto const p : id_val_pairs) {
		std::size_t const s(pnts->size());
		pnts->push_back(new GeoLib::Point(*(mesh.getNodes()[p.first]),
			mesh.getNodes()[p.first]->getID()));
		name_id_map->insert(std::make_pair("p"+std::to_string(p.first),s));
	}
	GeoLib::GEOObjects geo_obj;
	std::string geometry_name("SurfaceNodes");
	geo_obj.addPointVec(pnts, geometry_name, name_id_map);
	FileIO::BoostXmlGmlInterface io(geo_obj);
	io.setNameForExport(geometry_name);
	io.writeToFile(BaseLib::dropFileExtension(out_fname)+".gml");
}

void writeNeumannBCs(
	std::vector<std::pair<std::size_t,double>> const& id_val_pairs,
	std::string const& out_fname)
{
	std::ofstream st_out(BaseLib::dropFileExtension(out_fname)+".st");
	for (std::size_t k(0); k<id_val_pairs.size(); ++k) {
		st_out << "#SOURCE_TERM\n";
		st_out << "  $PCS_TYPE\n";
		st_out << "    GROUNDWATER_FLOW\n";
		st_out << "  $PRIMARY_VARIABLE\n";
		st_out << "    HEAD\n";
		st_out << "  $DIS_TYPE\n";
		st_out << "    CONSTANT " << id_val_pairs[k].second << "\n";
		st_out << "  $GEO_TYPE\n";
		st_out << "    POINT p" << id_val_pairs[k].first << "\n";
		st_out << "  $TIM_TYPE\n";
		st_out << "    CURVE " << k+1 << "\n";
	}
	st_out << "#STOP\n";
	st_out.close();
}

std::vector<std::pair<std::size_t, double>>
computeSurfaceNodeIDAndArea(MeshLib::Mesh const& mesh)
{
	const MathLib::Vector3 dir(0,0,-1);
	double const angle(80);
	std::unique_ptr<MeshLib::Mesh> sfc_mesh(
		MeshLib::MeshSurfaceExtraction::getMeshSurface(mesh, dir, angle, true));

	if (sfc_mesh.get() == nullptr) {
		WARN("Could not extract surface mesh.");
		return std::vector<std::pair<std::size_t,double>>();
	}

	std::vector<double> const area_vec =
		MeshLib::MeshSurfaceExtraction::getSurfaceAreaForNodes(*sfc_mesh);
	// put id info and area info into pair
	std::vector<std::pair<std::size_t, double>> node_area_vec;
	for (std::size_t k(0); k<area_vec.size(); ++k) {
		node_area_vec.push_back(
			std::make_pair(sfc_mesh->getNodes()[k]->getID(), area_vec[k]));
	}
	return node_area_vec;
}

std::vector<std::pair<size_t,double>>
computeNeumannConditionsWithSurfaceIntegration(
	MeshLib::Mesh const& subsurface_mesh,
	std::vector<std::pair<std::size_t, double>> const& node_area_vec,
	GeoLib::Raster const& raster, double scaling)
{
	double const no_data(raster.getNoDataValue());
	std::vector<std::pair<std::size_t,double>> direct_values;
	direct_values.reserve(node_area_vec.size());
	for (std::size_t i=0; i<node_area_vec.size(); ++i) {
		std::size_t const id(node_area_vec[i].first);
		double val(raster.getValueAtPoint(*(subsurface_mesh.getNodes()[id])));
		if (val == no_data)
			val = 0.0;
		direct_values.push_back(std::make_pair(id, val));
	}

	return direct_values;
}

int main(int argc, char *argv[])
{
	LOGOG_INITIALIZE();
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog::Cout *logogCout(new logog::Cout);
	logogCout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Creates OGS-5 transient Neumann conditions", ' ', BaseLib::BuildInfo::git_describe);
	TCLAP::ValueArg<std::string> input_arg("i", "input",
		"file containing a mesh (*.msh or *.vtu)", true, "", "string");
	TCLAP::UnlabeledMultiArg<std::string> raster_input_arg("input-raster",
		"raster file (*.asc) that will be used for the computation",
		true, "string");
	TCLAP::ValueArg<std::string> output_arg("o", "output-name",
		"filename for output", true, "", "string");
	cmd.add(input_arg);
	cmd.add(raster_input_arg);
	cmd.add(output_arg);
	cmd.parse(argc, argv);

	std::string const& mesh_fname(input_arg.getValue());
	std::string const& out_fname(output_arg.getValue());

	// load mesh from file
	std::unique_ptr<MeshLib::Mesh> mesh(FileIO::readMeshFromFile(mesh_fname));
	if (mesh.get() == nullptr)
		return 1;

	// preprocessing: compute surface node ids and the corresponding surface area
	std::vector<std::pair<std::size_t, double>> area_per_surface_node =
		computeSurfaceNodeIDAndArea(*(mesh.get()));

	// write geometry and the neumann file referencing the curves that will be
	// written below
	writeGeometry(*(mesh.get()), area_per_surface_node, out_fname);
	writeNeumannBCs(area_per_surface_node, out_fname);

	std::vector<std::vector<double>> neumann_vals(area_per_surface_node.size());
	// load and process raster files
	for (auto it=raster_input_arg.begin(); it != raster_input_arg.end(); ++it) {
		INFO("processing raster file \"%s\".", it->c_str());
		std::unique_ptr<GeoLib::Raster> raster(
			FileIO::AsciiRasterInterface::getRasterFromASCFile(*it));
			std::vector<std::pair<std::size_t, double>> result(
				computeNeumannConditionsWithSurfaceIntegration(
					*(mesh.get()), area_per_surface_node, *(raster.get()), 1.0));
		if (!result.empty()) {
			// transmit the results to the global data structure
			for (std::size_t k(0); k<result.size(); ++k) {
				neumann_vals[k].push_back(result[k].second);
			}
		}
	}

	boost::gregorian::date d0(1977, 12, 1);
	// write tim and neumann boundary condition to file
	std::ofstream out(BaseLib::dropFileExtension(out_fname)+".rfd");
	std::ofstream out_tim(BaseLib::dropFileExtension(out_fname)+".tim");
	out_tim << "#TIME_STEPPING\n";
	out_tim << "  $PCS_TYPE\n";
	out_tim << "    GROUNDWATER_FLOW\n";
	out_tim << "  $TIME_START\n";
	out_tim << "    0\n";
	out_tim << "  $TIME_END\n";
	out_tim << "    1e35\n";
	out_tim << "  $TIME_STEPS\n";
	for (std::size_t k(0); k<neumann_vals.size(); ++k) {
		out << "#CURVE ; " << area_per_surface_node[k].first << "\n";
		for (std::size_t j(0); j<neumann_vals[k].size(); ++j) {
			boost::gregorian::date dj(d0 + boost::gregorian::months(j));
			// for transfering values from meter per month to meter per second
			boost::gregorian::days days_per_month(d0+boost::gregorian::months(j+1)-dj);
			boost::gregorian::days days(dj-d0);
			out << days.days()*86400 << " "
				<< neumann_vals[k][j] / (days_per_month.days() * 86400 * 1000) << "\n";
			if (k == 0) {
				boost::gregorian::days days_diff(d0+boost::gregorian::months(j+1)-dj);
				out_tim << "1 " << days_diff.days() * 86400 << "\n";
			}
		}
	}
	out_tim << "#STOP\n";
	out << "#STOP\n";
	out_tim.close();
	out.close();

	delete logogCout;
	delete custom_format;
	LOGOG_SHUTDOWN();
	return 0;
}

