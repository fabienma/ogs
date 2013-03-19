/**
 * @file GocadSGridReader.cpp
 * @author Thomas Fischer
 * @date Mar 7, 2013
 * @brief 
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "GocadSGridReader.h"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>


// boost
#include <boost/tokenizer.hpp>

// ThirdParty/logog
#include "logog/include/logog.hpp"

// GeoLib
#include "GEOObjects.h"

// MeshLib
#include "Elements/Hex.h"

namespace FileIO
{

GocadSGridReader::GocadSGridReader(std::string const& fname) :
		_fname(fname), _path(_fname.substr(0, _fname.find_last_of("/\\") + 1))
{
	// check if file exists
	std::ifstream in(_fname.c_str());
	if (!in) {
		ERR("Could not open \"%s\".", _fname.c_str());
		in.close();
		return;
	}

	// read information in the stratigraphic grid file
	std::string line;
	while (std::getline(in, line)) {
		if (line.find("AXIS_N ") != std::string::npos) {
			parseDims(line);
		}
		if (line.find("POINTS_FILE ") != std::string::npos) {
			parsePointsFileName(line);
		}
		if (line.find("PROP_FILE ") != std::string::npos) {
			parsePropertiesFileName(line);
		}
//		if (line.find("FLAGS_FILE ") != std::string::npos
//				&& line.find("REGION_FLAGS_FILE ") == std::string::npos) {
//			parseFlagsFileName(line);
//		}
//		if (line.find("REGION_FLAGS_FILE ") != std::string::npos) {
//			parseRegionFlagsFileName(line);
//		}
	}
	in.close();

	readNodesBinary();
	makeNodesUnique();
	readElementPropertiesBinary();
	createElements();
	readSplitNodesAndModifyElements();
	removeNullVolumeElements();
}

GocadSGridReader::~GocadSGridReader()
{
}

void GocadSGridReader::parseDims(std::string const& line)
{
	std::size_t x_dim(0), y_dim(0), z_dim(0);
	boost::tokenizer<> tok(line);
	auto it(tok.begin());
	it++; // overread token "AXIS"
	it++; // overread "N"
	std::stringstream ssx(*(it), std::stringstream::in | std::stringstream::out);
	ssx >> x_dim;
	it++;
	std::stringstream ssy(*it, std::stringstream::in | std::stringstream::out);
	ssy >> y_dim;
	it++;
	std::stringstream ssz(*it, std::stringstream::in | std::stringstream::out);
	ssz >> z_dim;
	_index_calculator = IndexCalculator(x_dim, y_dim, z_dim);
}

void GocadSGridReader::parsePointsFileName(std::string const& line)
{
	std::size_t beg_pos(line.find_first_of(" ") + 1);
	_pnts_fname = _path + line.substr(beg_pos, line.length() - beg_pos);
}

void GocadSGridReader::parsePropertiesFileName(std::string const& line)
{
	std::size_t beg_pos(line.find_first_of(" ") + 3);
	_properties_fname = _path + line.substr(beg_pos, line.length() - beg_pos);
}

void GocadSGridReader::parseFlagsFileName(std::string const& line)
{
	std::size_t beg_pos(line.find_first_of(" ") + 1);
	_flags_fname = _path + line.substr(beg_pos, line.length() - beg_pos);
}

void GocadSGridReader::parseRegionFlagsFileName(std::string const& line)
{
	std::size_t beg_pos(line.find_first_of(" ") + 1);
	_region_flags_fname = _path + line.substr(beg_pos, line.length() - beg_pos);
}

void GocadSGridReader::readNodesBinary()
{
	std::ifstream in(_pnts_fname.c_str(), std::ios::in | std::ios::binary);
	if (!in) {
		ERR("Could not open points file \"%s\".", _pnts_fname.c_str());
		return;
	}

	_nodes.resize(_index_calculator._n_nodes);

	char inbuff[12], reword[4];
	double coords[3];

	for (std::size_t k(0); k < _index_calculator._n_nodes; k++) {
		in.read((char*) inbuff, 12 * sizeof(char));
		for (std::size_t i = 0; i < 12; i += 4) {
			for (std::size_t j = 0; j < 4; j++) {
				reword[j] = inbuff[i + 3 - j];
			}
			coords[i/4] = static_cast<double>(*reinterpret_cast<float*>(&reword[0]));
		}
		_nodes[k] = new MeshLib::Node(coords, k);
	}
	in.close();
}

void GocadSGridReader::readElementPropertiesBinary()
{
	std::ifstream in(_properties_fname.c_str());
	if (!in) {
		ERR("Could not open element property file \"%s\".", _properties_fname.c_str());
		in.close();
		return;
	}

	_properties.resize(_index_calculator._n_cells);
	char inbuff[4], reword[4];

	const std::size_t n_properties(_properties.size());
	for (std::size_t k(0); k < n_properties; k++) {
		in.read((char*) inbuff, 4 * sizeof(char));
		if (in) {
			for (std::size_t j = 0; j < 4; j++) {
				reword[j] = inbuff[3 - j];
			}
			_properties[k] = static_cast<double>(*reinterpret_cast<float*>(&reword[0]));
		} else {
			k = n_properties;
		}
	}
}

//void GocadSGridReader::readFlagsBinary(std::size_t n_flags,
//		std::vector<double> &flags)
//{
//	std::ifstream in(_flags_fname.c_str());
//	if (!in) {
//		std::cout << "Could not open " << _flags_fname << "." << std::endl;
//		in.close();
//		return;
//	}
//
//	char inbuff[4], reword[4];
//
//	for (std::size_t k(0); k < n_flags; k++) {
//		in.read((char*) inbuff, 4 * sizeof(char));
//		for (std::size_t i = 0; i < 4; i += 4) {
//			for (std::size_t j = 0; j < 4; j++) {
//				reword[j] = inbuff[i + 3 - j];
//			}
//			if (k > n_flags - 20)
//				std::cout << (*reinterpret_cast<int*>(&reword[0])) << std::endl;
//			flags[4 * k + i / 4] = static_cast<double>(*reinterpret_cast<float*>(&reword[0]));
//		}
//	}
//}
//
//void readBinarySGridRegionFlags(std::string const& fname, std::size_t n_region_flags,
//		std::vector<double> &region_flags)
//{
//	std::ifstream in(fname.c_str());
//	if (!in) {
//		std::cout << "Could not open " << fname << "." << std::endl;
//		in.close();
//		return;
//	}
//
//	char inbuff[16], reword[4];
//
//	for (std::size_t k(0); k < n_region_flags; k++) {
//		in.read((char*) inbuff, 16 * sizeof(char));
//		for (std::size_t i = 0; i < 16; i += 4) {
//			for (std::size_t j = 0; j < 4; j++) {
//				reword[j] = inbuff[i + 3 - j];
//			}
//			region_flags[4 * k + i / 4] =
//					static_cast<double>(*reinterpret_cast<float*>(&reword[0]));
//		}
//	}
//}

void GocadSGridReader::makeNodesUnique()
{
	// to make nodes unique GeoLib::PointVec will be employed
	std::string pname("test");
	std::vector<GeoLib::Point*>* tmp_nodes(new std::vector<GeoLib::Point*>(_nodes.size()));
	std::copy(_nodes.begin(), _nodes.end(), tmp_nodes->begin());
	_nodes.clear();

	GeoLib::GEOObjects geo;
	geo.addPointVec(tmp_nodes, pname);

	// save node <-> id mapping needed for creating the mesh elements later on
	_node_id_map.resize(_index_calculator._n_nodes);
	std::vector<std::size_t> const& node_id_map(geo.getPointVecObj(pname)->getIDMap());
	std::copy(node_id_map.begin(), node_id_map.end(), _node_id_map.begin());

	// copy unique nodes
	tmp_nodes = const_cast<std::vector<GeoLib::Point*>*>(geo.getPointVecObj(pname)->getVector());
	_nodes.resize(tmp_nodes->size());
	for (std::size_t k(0); k < tmp_nodes->size(); k++) {
		_nodes[k] = new MeshLib::Node((*tmp_nodes)[k]->getCoords());
	}
}

void GocadSGridReader::createElements()
{
	_elements.resize(_index_calculator._n_cells);
	std::array<MeshLib::Node*, 8> element_nodes;
	std::size_t cnt(0);
	for (std::size_t k(0); k < _index_calculator._z_dim-1; k++) {
		for (std::size_t j(0); j < _index_calculator._y_dim-1; j++) {
			for (std::size_t i(0); i < _index_calculator._x_dim-1; i++) {
				element_nodes[0] = _nodes[_node_id_map[_index_calculator(i,j,k)]];
				element_nodes[1] = _nodes[_node_id_map[_index_calculator(i+1,j,k)]];
				element_nodes[2] = _nodes[_node_id_map[_index_calculator(i+1,j+1,k)]];
				element_nodes[3] = _nodes[_node_id_map[_index_calculator(i,j+1,k)]];
				element_nodes[4] = _nodes[_node_id_map[_index_calculator(i,j,k+1)]];
				element_nodes[5] = _nodes[_node_id_map[_index_calculator(i+1,j,k+1)]];
				element_nodes[6] = _nodes[_node_id_map[_index_calculator(i+1,j+1,k+1)]];
				element_nodes[7] = _nodes[_node_id_map[_index_calculator(i,j+1,k+1)]];
				_elements[cnt] = new MeshLib::Hex(element_nodes, static_cast<unsigned>(_properties[_index_calculator(i,j,k)]));
				cnt++;
			}
		}
	}
}

void GocadSGridReader::readSplitNodesAndModifyElements()
{
	std::ifstream in(_fname.c_str());
	if (!in) {
		ERR("Could not open \"%s\".", _fname.c_str());
		in.close();
		return;
	}

	// read information in the stratigraphic grid file
	std::string line;
	std::stringstream ss;
	while (std::getline(in, line)) {
		std::size_t pos(line.find("SPLIT "));
		if (pos != std::string::npos) {
			ss << line.substr(pos+6, line.size()-(pos+6));
			// read position in grid
			std::size_t u, v, w;
			ss >> u;
			ss >> v;
			ss >> w;
			// read coordinates for the split node
			double coords[3];
			ss >> coords[0];
			ss >> coords[1];
			ss >> coords[2];
			// read the id
			std::size_t id;
			ss >> id;
			// read the affected cells
			std::array<bool, 8> cells;
			for (std::size_t k(0); k<cells.size(); k++) {
				ss >> cells[k];
			}
			std::size_t new_node_pos(_nodes.size());
			_nodes.push_back(new MeshLib::Node(coords));
			// get mesh node to substitute in elements
			MeshLib::Node const*const node2sub(_nodes[_node_id_map[_index_calculator(u,v,w)]]);
			if (cells[0]) {
				modifyElement(u, v, w, node2sub, _nodes[new_node_pos]);
			}
			if (cells[1]) {
				modifyElement(u - 1, v, w, node2sub, _nodes[new_node_pos]);
			}
			if (cells[2]) {
				modifyElement(u, v - 1, w, node2sub, _nodes[new_node_pos]);
			}
			if (cells[3]) {
				modifyElement(u - 1, v - 1, w, node2sub, _nodes[new_node_pos]);
			}
			if (cells[4]) {
				modifyElement(u, v, w - 1, node2sub, _nodes[new_node_pos]);
			}
			if (cells[5]) {
				modifyElement(u - 1, v, w - 1, node2sub, _nodes[new_node_pos]);
			}
			if (cells[6]) {
				modifyElement(u, v - 1, w - 1, node2sub, _nodes[new_node_pos]);
			}
			if (cells[7]) {
				modifyElement(u - 1, v - 1, w - 1, node2sub, _nodes[new_node_pos]);
			}
		}
	}
}

void GocadSGridReader::modifyElement(std::size_t u, std::size_t v, std::size_t w,
		MeshLib::Node const* node2sub, MeshLib::Node * substitute_node)
{
	// ensure (u,v,w) is a valid cell
	if (u >= _index_calculator._x_dim - 1 || v >= _index_calculator._y_dim - 1
			|| w >= _index_calculator._z_dim - 1) {
		return;
	}

	// get the hex in which the node should be substituted
	const std::size_t idx(_index_calculator.getCellIdx(u,v,w));
	if (idx == std::numeric_limits<std::size_t>::max())
		return;

	MeshLib::Element *hex(_elements[idx]);
	// get the node pointers of the cell
	MeshLib::Node *const* hex_nodes(hex->getNodes());
	// search for the position the split node will be set to
	MeshLib::Node *const* node_pos(std::find(hex_nodes, hex_nodes+8, node2sub));
	// set the split node instead of the node2sub
	if (node_pos != hex_nodes+8) {
		const_cast<MeshLib::Node**>(hex_nodes)[std::distance(hex_nodes, node_pos)] = substitute_node;
	}
}

void GocadSGridReader::removeNullVolumeElements()
{
	auto const new_end(
			std::remove_if(_elements.begin(), _elements.end(),
					[](MeshLib::Element *elem) {
						return (elem->getContent() < std::numeric_limits<double>::epsilon());
					}
			)
	);
	_elements.erase(new_end, _elements.end());
}

} // end namespace FileIO
