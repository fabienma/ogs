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

// BaseLib
#include "FileTools.h"

// GeoLib
#include "AABB.h"

// MeshLib
#include "Elements/Hex.h"
#include "Elements/Quad.h"

namespace FileIO
{

typedef boost::dynamic_bitset<> Bitset;

typedef GocadSGridReader::Region Region;
typedef GocadSGridReader::Layer Layer;
typedef GocadSGridReader::GocadProperty GocadProperty;

std::ostream& operator<<(std::ostream& os, Region const& r)
{
	return os << "(" << r.name << "|" << r.bit << ")";
}

std::ostream& operator<<(std::ostream& os, Layer const& l)
{
	std::copy(l.regions.begin(), l.regions.end(),
		std::ostream_iterator<Region>(os, " "));
	return os;
}

std::ostream& operator<<(std::ostream& os, GocadProperty const& p)
{
	return os << p._property_name << " " << p._property_id << " "
			<< p._property_data_type << " " << p._property_data_fname;
}

Region parseRegion(std::string const& line)
{
	std::istringstream iss(line);
	std::istream_iterator<std::string> it(iss);
	// Check first word is REGION or MODEL_REGION.
	if (*it != std::string("REGION") && *it != std::string("MODEL_REGION"))
	{
		ERR("Expected REGION or MODEL_REGION keyword but \"%s\" found.\n", it->c_str());
		throw std::runtime_error("In parseRegion() expected REGION or MODEL_REGION keyword not found.\n");
	}
	++it;

	Region r;
	r.name = *it;
	++it;
	r.bit = atoi(it->c_str());

	return r;
}

Layer parseLayer(std::string const& line, std::vector<Region> const& regions)
{
	std::istringstream iss(line);
	std::istream_iterator<std::string> it(iss);
	// Check first word is MODEL_LAYER.
	if (*it != std::string("MODEL_LAYER"))
	{
		ERR("Expected MODEL_LAYER keyword but \"%s\" found.\n", it->c_str());
		throw std::runtime_error("In parseRegion() expected MODEL_LAYER keyword not found.\n");
	}
	++it;

	Layer l;
	while (it != std::istream_iterator<std::string>() && *it != "END")
	{
		l.regions.push_back(
			*std::find_if(regions.begin(), regions.end(),
				[&](Region const& r) { return r.name == *it; }));
		++it;
	}

	return l;
}

GocadProperty parseGocadPropertyMetaData(std::string &line, std::istream &in, std::string const& path)
{
	boost::char_separator<char> sep("\t ");
	boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
	auto tok_it(tokens.begin());
	// A property section in Gocad file starts with a line
	// PROPERTY id "property_name"
	if (tok_it->compare("PROPERTY") != 0 ) {
		ERR("Expected PROPERTY keyword but \"%s\" found.", tok_it->c_str());
		throw std::runtime_error("In parseGocadPropertyMetaData() expected PROPERTY keyword not found.");
	}
	tok_it++;

	GocadProperty prop;
	prop._property_id = BaseLib::str2number<std::size_t>(*tok_it);
	tok_it++;
	prop._property_name = *tok_it;
	BaseLib::trim(prop._property_name, '\"');

	while (getline(in, line)) {
		tokens.assign(line);

		tok_it = tokens.begin();
		// this is the last entry of the property
		if (tok_it->compare("PROP_FILE") == 0) {
			tok_it++;
			if (! prop.checkID(*tok_it)) {
				throw std::runtime_error("parseGocadPropertyMetaData(): id mismatch.");
			}
			tok_it++;
			std::string tmp(*tok_it);
			BaseLib::trim(tmp, '\"');
			prop._property_data_fname = path + tmp;
			return prop;
		}

		if (tok_it->compare("PROPERTY_CLASS") == 0) {
			tok_it++;
			if (! prop.checkID(*tok_it)) {
				throw std::runtime_error("parseGocadPropertyMetaData(): id mismatch.");
			}
			tok_it++;
			prop._property_class_name = *tok_it;
		}

		if (tok_it->compare("PROPERTY_CLASS") == 0) {
			tok_it++;
			if (! prop.checkID(*tok_it)) {
				throw std::runtime_error("parseGocadPropertyMetaData(): id mismatch.");
			}
			tok_it++;
			prop._property_class_name = *tok_it;
		}

		if (tok_it->compare("PROPERTY_SUBCLASS") == 0) {
			tok_it++;
			if (! prop.checkID(*tok_it)) {
				throw std::runtime_error("parseGocadPropertyMetaData(): id mismatch.");
			}
			tok_it++;
			if (tok_it->compare("QUANTITY") != 0) {
				ERR("Expected keyword QUANTITY but found \"%s\".", tok_it->c_str());
				throw std::runtime_error("parseGocadPropertyMetaData(): id mismatch.");
			}
			tok_it++;
			prop._property_data_type = *tok_it;
		}

		if (tok_it->compare("PROP_UNIT") == 0 ||
				tok_it->compare("PROP_ORIGINAL_UNIT") == 0) {
			tok_it++;
			if (! prop.checkID(*tok_it)) {
				throw std::runtime_error("parseGocadPropertyMetaData(): id mismatch.");
			}
			tok_it++;
			prop._property_unit = *tok_it;
		}

		if (tok_it->compare("PROP_NO_DATA_VALUE") == 0) {
			tok_it++;
			if (! prop.checkID(*tok_it)) {
				throw std::runtime_error("parseGocadPropertyMetaData(): id mismatch.");
			}
			tok_it++;
			prop._property_no_data_value = BaseLib::str2number<double>(*tok_it);
		}
	}
	return prop;
}

GocadSGridReader::GocadSGridReader(std::string const& fname) :
		_fname(fname), _path(_fname.substr(0, _fname.find_last_of("/\\") + 1)),
		_n_face_sets(0), _double_precision_binary(false), _bin_pnts_in_double_precision(false)

{
	// check if file exists
	std::ifstream in(_fname.c_str());
	if (!in) {
		ERR("Could not open \"%s\".", _fname.c_str());
		in.close();
		return;
	}

	bool pnts_read(false);

	// read information in the stratigraphic grid file
	std::string line;
	while (std::getline(in, line)) {
		if (line.compare(0, 8, "HEADER {") == 0)
		{
			parseHeader(in);
		}
		if (line.compare(0, 7, "AXIS_N ") == 0)
		{
			parseDims(line);
		}
		else if (line.compare(0, 12, "POINTS_FILE ") == 0)
		{
			parseFileName(line, _pnts_fname);
		}
		else if (line.compare(0, 9, "PROPERTY ") == 0)
		{
			_property_meta_data_vecs.push_back(parseGocadPropertyMetaData(line, in, _path));
		}
		else if (line.compare(0, 35, "BINARY_POINTS_IN_DOUBLE_PRECISION 1") == 0) {
			_bin_pnts_in_double_precision = true;
		}
		else if (line.compare(0, 11, "FLAGS_FILE ") == 0) {
			parseFileName(line, _flags_fname);
		}
		else if (line.compare(0, 18, "REGION_FLAGS_FILE ") == 0)
		{
			parseFileName(line, _region_flags_fname);
		}
		else if (line.compare(0, 7, "REGION ") == 0 || line.compare(0, 13, "MODEL_REGION ") == 0)
		{
			regions.push_back(parseRegion(line));
		}
		else if (line.compare(0, 12, "MODEL_LAYER ") == 0)
		{
			layers.push_back(parseLayer(line, regions));
		}
		else if (line.compare(0, 24, "REGION_FLAGS_BIT_LENGTH ") == 0)
		{
			std::istringstream iss(line);
			std::istream_iterator<std::string> it(iss);
			it++;
			std::size_t bit_length = atoi(it->c_str());
			if (regions.size() != bit_length)
			{
				ERR("%d regions read but %d expected.\n", regions.size(), bit_length);
				throw std::runtime_error("Number of read regions differs from expected.\n");
			}
		}
		else if (line.compare(0, 9, "FACE_SET ") == 0) {
			// first read the points
			if (! pnts_read) {
				readNodesBinary();
				pnts_read = true;
			}
			parseFaceSet(line, in);
		}
		else
		{
			//std::cout << "Skip: \"" << line << "\"\n";
		}
	}
	std::cout << regions.size() << " regions read:\n";
	std::copy(regions.begin(), regions.end(), std::ostream_iterator<Region>(std::cout, "\t"));
	std::cout << "\n";
	std::cout << layers.size() << " layers read:\n";
	std::copy(layers.begin(), layers.end(), std::ostream_iterator<Layer>(std::cout, "\n"));

	std::cout << "meta data for " << _property_meta_data_vecs.size() << " properties read:\n";
	std::copy(_property_meta_data_vecs.begin(), _property_meta_data_vecs.end(),
			std::ostream_iterator<GocadProperty>(std::cout, "\n"));

	// if not done already read the points
	if (! pnts_read) {
		readNodesBinary();
		pnts_read = true;
	}
	readElementPropertiesBinary();
	std::vector<Bitset> region_flags = readRegionFlagsBinary();
	//mapRegionFlagsToCellProperties(region_flags);	// modifies _material_ids.

	createElements();
	readSplitNodesAndModifyElements();

	GeoLib::AABB<MeshLib::Node> aabb(_nodes.begin(), _nodes.end());
	MeshLib::Node center_node((aabb.getMaxPoint()[0] + aabb.getMinPoint()[0])/2.0,
			(aabb.getMaxPoint()[1] + aabb.getMinPoint()[1])/2.0, 0.0);
	INFO("translated model (-%f, -%f, -%f).", center_node[0], center_node[1], center_node[2]);
	std::for_each(_nodes.begin(), _nodes.end(),
			[&center_node](MeshLib::Node* node)
			{
				(*node)[0] -= center_node[0];
				(*node)[1] -= center_node[1];
			}
	);

//	removeNullVolumeElements();
	GocadProperty face_set_property;
	face_set_property._property_id = 0;
	face_set_property._property_name = "CellIDs";
	face_set_property._property_class_name = "CellIDsData";
	face_set_property._property_unit = "unitless";
	face_set_property._property_data_type = "double";
	face_set_property._property_data_fname = "";
	face_set_property._property_no_data_value = -1.0;
	face_set_property._property_data.resize(_index_calculator._n_cells);
	std::iota(face_set_property._property_data.begin(), face_set_property._property_data.end(),
			0.0);

	_property_meta_data_vecs.push_back(face_set_property);
	in.close();
}

GocadSGridReader::~GocadSGridReader()
{
}

void GocadSGridReader::parseHeader(std::istream &in)
{
	std::string line;
	while (std::getline(in, line)) {
		if (line.compare(0, 1, "}") == 0) {
			return;
		}
		if (line.compare(0, 27, "double_precision_binary: on") == 0) {
			_double_precision_binary = true;
		}
	}
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

void GocadSGridReader::parseFileName(std::string const& line, std::string &result_string) const
{
	std::size_t beg_pos(line.find_first_of(" ") + 1);
	std::string fname(line.substr(beg_pos, line.length() - beg_pos));
	BaseLib::trim(fname, '\"');
	result_string = _path + fname;
}

/**
 * @param line input/output
 * @param in input stream containing the face set
 */
void GocadSGridReader::parseFaceSet(std::string &line, std::istream &in)
{
	// create and initialize a GocadProperty object for storing face set data
	GocadProperty face_set_property;
	face_set_property._property_id = _n_face_sets;
	face_set_property._property_name = "FaceSet";
	face_set_property._property_class_name = "FaceSetData";
	face_set_property._property_unit = "unitless";
	face_set_property._property_data_type = "double";
	face_set_property._property_data_fname = "";
	face_set_property._property_no_data_value = -1.0;
	face_set_property._property_data.resize(_index_calculator._n_cells);
	std::fill(face_set_property._property_data.begin(), face_set_property._property_data.end(),
			face_set_property._property_no_data_value);

	std::istringstream iss(line);
	std::istream_iterator<std::string> it(iss);
	// Check first word is FACE_SET
	if (*it != std::string("FACE_SET")) {
		ERR("Expected FACE_SET keyword but \"%s\" found.", it->c_str());
		throw std::runtime_error("In GocadSGridReader::parseFaceSet() expected FACE_SET keyword not found.");
	}
	++it;
	face_set_property._property_name += BaseLib::number2str(*it);
	++it;
	std::size_t const n_of_face_set_ids(static_cast<std::size_t>(atoi(it->c_str())));
	std::size_t face_set_id_cnt(0);

	while (getline(in, line) && face_set_id_cnt < n_of_face_set_ids) {
		boost::char_separator<char> sep("\t ");
		boost::tokenizer<boost::char_separator<char> > tokens(line, sep);

		for(auto tok_it  = tokens.begin(); tok_it != tokens.end(); ) {
			std::size_t id(static_cast<std::size_t>(atoi(tok_it->c_str())));
			tok_it++;
			std::size_t face_indicator(static_cast<std::size_t>(atoi(tok_it->c_str())));
			tok_it++;

			if (id >= _index_calculator._n_nodes) {
				ERR("Face set id %d is greater than the number of nodes (%d).", id, _index_calculator._n_nodes);
			} else {
				dynamic_cast<MeshLib::GocadNode*>(_nodes[id])->setFaceSet(_n_face_sets+1, face_indicator);
				std::array<std::size_t,3> coords(_index_calculator.getCoordsForID(id));
				if (coords[0] >= _index_calculator._x_dim-1)
					ERR("****** i coord %d to big for id %d.", coords[0], id);
				if (coords[1] >= _index_calculator._y_dim-1)
					ERR("****** j coord %d to big for id %d.", coords[1], id);
				if (coords[2] >= _index_calculator._z_dim-1)
					ERR("****** k coord %d to big for id %d.", coords[2], id);
				std::size_t const cell_id( _index_calculator.getCellIdx(coords[0], coords[1], coords[2]));
				face_set_property._property_data[cell_id] = face_indicator;
			}
			face_set_id_cnt++;
		}
	}

	if (face_set_id_cnt != n_of_face_set_ids) {
		ERR("Expected %d number of face set ids, read %d.", n_of_face_set_ids, face_set_id_cnt);
		throw std::runtime_error("Expected number of face set points does not match number of read points.");
	}
	_n_face_sets++;

	_property_meta_data_vecs.push_back(face_set_property);
}

// Reads given number of bits (rounded up to next byte) into a bitset.
// Used for reading region information which can be represented by some
// number of bits.
Bitset readBits(std::ifstream& in, const std::size_t bits)
{
	typedef Bitset::block_type block_t;
	std::size_t const bytes = static_cast<std::size_t>(std::ceil(bits/8.));
	std::size_t const blocks = (bytes + 1)/ sizeof(block_t);

	block_t data[blocks];
	std::fill_n(data, blocks, 0);
	in.read(reinterpret_cast<char*>(data), bytes);

	return Bitset(data, data + blocks);
}

void GocadSGridReader::readNodesBinary()
{
	std::ifstream in(_pnts_fname.c_str(), std::ios::binary);
	if (!in) {
		ERR("Could not open points file \"%s\".", _pnts_fname.c_str());
		throw std::runtime_error("Could not open points file.");
	}

	std::size_t const n = _index_calculator._n_nodes;
	_nodes.resize(n);

	double coords[3];

	std::size_t k = 0;
	while (in && k < n * 3)
	{
		if (_bin_pnts_in_double_precision) {
			coords[k % 3] = BaseLib::swapEndianness(BaseLib::readBinaryValue<double>(in));
		} else {
			coords[k % 3] = BaseLib::swapEndianness(BaseLib::readBinaryValue<float>(in));
		}
		if ((k + 1) % 3 == 0)
			_nodes[k/3] = new MeshLib::GocadNode(coords, k/3);
		k++;
	}
	if (k != n * 3 && !in.eof())
		ERR("Read different number of points. Expected %d floats, got %d.\n", n * 3, k);
}

void GocadSGridReader::mapRegionFlagsToCellProperties(std::vector<Bitset> const& rf)
{
//	std::size_t const n = _index_calculator._n_cells;
//	_material_ids.resize(n);
//	std::fill(_material_ids.begin(), _material_ids.end(), -1);
//	// region flags are stored in each node ijk and give the region index for the
//	// ijk-th cell.
//	for (std::size_t k(0); k < _index_calculator._z_dim-1; k++) {
//		for (std::size_t j(0); j < _index_calculator._y_dim-1; j++) {
//			for (std::size_t i(0); i < _index_calculator._x_dim-1; i++) {
//				// Find layers containing regions given by bits.
//				// Run over bits and push back set bits
//				std::set<std::size_t> layers_set;
//				for (auto r = regions.begin(); r != regions.end(); ++r)
//				{
//					if (rf[_index_calculator(i,j,k)].test(r->bit))
//					{
//						// Bit is set, find a layer.
//						for (std::size_t l_id = 0; l_id < layers.size(); ++l_id)
//							if (layers[l_id].hasRegion(*r))
//								layers_set.insert(l_id);
//					}
//				}
//				if (layers_set.size() != 1)
//					ERR("A cell %d %d %d belongs to multiple (%d) layers.", i, j, k, layers_set.size());
//
//				_material_ids[_index_calculator.getCellIdx(i,j,k)] = *layers_set.begin();
//			}
//		}
//	}
}

void GocadSGridReader::readElementPropertiesBinary()
{
	for (auto prop_it(_property_meta_data_vecs.begin());
			prop_it != _property_meta_data_vecs.end();
			prop_it++) {
		std::string const& fname(prop_it->_property_data_fname);
		if (prop_it->_property_data_fname.empty()) {
			WARN("Empty filename for property %s.", prop_it->_property_name.c_str());
			continue;
		}
		std::vector<float> float_properties =
				BaseLib::readBinaryArray<float>(fname, _index_calculator._n_cells);
		std::transform(float_properties.cbegin(), float_properties.cend(), float_properties.begin(),
				[](float const& val) {
					return BaseLib::swapEndianness(val);
				});

		std::vector<double> properties;
		prop_it->_property_data.resize(float_properties.size());
		std::copy(float_properties.begin(), float_properties.end(), prop_it->_property_data.begin());
		if (prop_it->_property_data.empty()) {
			ERR("Reading of element properties file \"%s\" failed.", fname.c_str());
		}
	}
}

std::vector<int> GocadSGridReader::readFlagsBinary() const
{
	std::vector<int> result;
	if (!_double_precision_binary) {
		result = BaseLib::readBinaryArray<int32_t>(_flags_fname, _index_calculator._n_nodes);
		std::for_each(result.begin(), result.end(),
						[](int32_t& val) {
							BaseLib::swapEndianness(val);
						});
	} else {
		result = BaseLib::readBinaryArray<int>(_flags_fname, _index_calculator._n_nodes);
		std::for_each(result.begin(), result.end(),
						[](int& val) {
							BaseLib::swapEndianness(val);
						});
	}

	if (result.empty())
		ERR("Reading of flags file \"%s\" failed.", _flags_fname.c_str());

	return result;
}

std::vector<Bitset> GocadSGridReader::readRegionFlagsBinary() const
{
	std::vector<Bitset> result;

	std::ifstream in(_region_flags_fname.c_str());
	if (!in) {
		ERR("readRegionFlagsBinary(): Could not open file \"%s\" for input.\n", _region_flags_fname.c_str());
		in.close();
		return result;
	}

	std::size_t const n = _index_calculator._n_nodes;
	result.resize(n);

	std::size_t k = 0;
	while (in && k < n)
	{
		result[k++] = readBits(in, regions.size());
	}
	if (k != n && !in.eof())
		ERR("Read different number of values. Expected %d, got %d.\n", n, k);

	return result;
}

void GocadSGridReader::createElements()
{
	_elements.resize(_index_calculator._n_cells);
	std::array<MeshLib::Node*, 8> element_nodes;
	std::size_t cnt(0);
	for (std::size_t k(0); k < _index_calculator._z_dim-1; k++) {
		for (std::size_t j(0); j < _index_calculator._y_dim-1; j++) {
			for (std::size_t i(0); i < _index_calculator._x_dim-1; i++) {
				element_nodes[0] = _nodes[_index_calculator(i,j,k)];
				element_nodes[1] = _nodes[_index_calculator(i+1,j,k)];
				element_nodes[2] = _nodes[_index_calculator(i+1,j+1,k)];
				element_nodes[3] = _nodes[_index_calculator(i,j+1,k)];
				element_nodes[4] = _nodes[_index_calculator(i,j,k+1)];
				element_nodes[5] = _nodes[_index_calculator(i+1,j,k+1)];
				element_nodes[6] = _nodes[_index_calculator(i+1,j+1,k+1)];
				element_nodes[7] = _nodes[_index_calculator(i,j+1,k+1)];
				_elements[cnt] = new MeshLib::Hex(element_nodes, _index_calculator.getCellIdx(i,j,k));
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

			static_cast<MeshLib::GocadNode*>(_nodes[_index_calculator(u,v,w)])->setSplit(true);
			std::size_t const new_node_pos(_nodes.size());
			MeshLib::GocadNode *new_node(new MeshLib::GocadNode(* static_cast<MeshLib::GocadNode*>(_nodes[_index_calculator(u,v,w)])));
			new_node->resetID(new_node_pos);
			(*new_node)[0] = coords[0];
			(*new_node)[1] = coords[1];
			(*new_node)[2] = coords[2];
			_nodes.push_back(new_node);

			// get mesh node to substitute in elements
			MeshLib::Node const*const node2sub(_nodes[_index_calculator(u,v,w)]);

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
						if (elem->getContent() < std::numeric_limits<double>::epsilon()) {
							delete elem;
							return true;
						}
						return false;
					}
			)
	);
	_elements.erase(new_end, _elements.end());

	std::size_t const n_elements(_elements.size());
	for (std::size_t k(0); k < n_elements; k++) {
		_elements[k]->setValue(k);
	}
}

boost::optional<GocadProperty const&>
GocadSGridReader::getProperty(std::string const& name) const
{
	auto const it(std::find_if(_property_meta_data_vecs.begin(), _property_meta_data_vecs.end(),
			[&name](GocadProperty const& p) { return p._property_name.compare(name) == 0; }));
	if (it != _property_meta_data_vecs.end())
		return boost::optional<GocadProperty const&>(*it);
	else
		return boost::optional<GocadProperty const&>();
}

std::vector<std::string> GocadSGridReader::getPropertyNames() const
{
	std::vector<std::string> names;
	std::transform(_property_meta_data_vecs.begin(),
			_property_meta_data_vecs.end(),
			std::back_inserter(names),
			[](GocadProperty const& p) {
				return p._property_name;
			}
		);
	return names;
}

MeshLib::Mesh* GocadSGridReader::getFaceSetMesh(std::size_t face_set_number) const
{
	auto const it = std::find_if(_property_meta_data_vecs.cbegin(), _property_meta_data_vecs.cend(),
			[face_set_number](GocadProperty const& p)
			{
				return (p._property_name.compare("FaceSet") == 0) && (p._property_id == face_set_number);
			}
		);
	if (it == _property_meta_data_vecs.cend())
		return nullptr;

	return nullptr;
}

} // end namespace FileIO
