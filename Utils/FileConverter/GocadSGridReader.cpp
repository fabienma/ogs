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

// MeshLib
#include "Elements/Hex.h"

namespace FileIO
{

typedef boost::dynamic_bitset<> Bitset;

typedef GocadSGridReader::Region Region;
typedef GocadSGridReader::Layer Layer;
typedef GocadSGridReader::FaceSet FaceSet;

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

/**
 *
 * @param line input/output
 * @param in input stream containing the face set
 * @return FaceSet
 */
FaceSet parseFaceSet(std::string &line, std::istream &in, std::size_t nu, std::size_t nv)
{
	std::istringstream iss(line);
	std::istream_iterator<std::string> it(iss);
	// Check first word is FACE_SET
	if (*it != std::string("FACE_SET")) {
		ERR("Expected FACE_SET keyword but \"%s\" found.", it->c_str());
		throw std::runtime_error("In parseFaceSet() expected FACE_SET keyword not found.");
	}
	++it;

	FaceSet fs;
	fs._name = *it;
	it++;
	std::size_t number_of_faces(static_cast<std::size_t>(atoi(it->c_str())));
	std::size_t faces_cnt(0);

	while (getline(in, line) && faces_cnt < number_of_faces) {
		boost::char_separator<char> sep("\t ");
		boost::tokenizer<boost::char_separator<char> > tokens(line, sep);

		for(auto tok_it  = tokens.begin(); tok_it != tokens.end(); ) {
			const std::size_t cell_id(static_cast<std::size_t>(atoi(tok_it->c_str())));
			const std::size_t u(cell_id/(nu*nv));
			const std::size_t v((cell_id%(nu*nv))/nu);
			const std::size_t w((cell_id%(nu*nv))%nu);
			tok_it++;
			const std::size_t face_dir(static_cast<std::size_t>(atoi(tok_it->c_str())));
			fs._face_pos_and_dir.push_back({u, v, w, face_dir});
			tok_it++;
			faces_cnt++;
		}
	}

	if (faces_cnt != number_of_faces) {
		ERR("Expected %d number of faces, read %d.", number_of_faces, faces_cnt);
		throw std::runtime_error("Expected number of faces does not match number of read faces.");
	}

	return fs;
}

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
		if (line.compare(0, 7, "AXIS_N ") == 0)
		{
			parseDims(line);
		}
		else if (line.compare(0, 12, "POINTS_FILE ") == 0)
		{
			parsePointsFileName(line);
		}
		else if (line.compare(0, 10, "PROP_FILE ") == 0)
		{
			parsePropertiesFileName(line);
		}
		else if (line.compare(0, 11, "FLAGS_FILE ") == 0) {
			parseFlagsFileName(line);
		}
		else if (line.compare(0, 18, "REGION_FLAGS_FILE ") == 0)
		{
			parseRegionFlagsFileName(line);
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
			_face_sets.push_back(parseFaceSet(line, in,
					_index_calculator._x_dim-1, _index_calculator._y_dim-1));
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

	readNodesBinary();

	readElementPropertiesBinary();
	std::vector<Bitset> region_flags = readRegionFlagsBinary();
	//mapRegionFlagsToCellProperties(region_flags);	// modifies _properties.

	createElements();
	readSplitNodesAndModifyElements();
	removeNullVolumeElements();

	in.close();
}

GocadSGridReader::~GocadSGridReader()
{
}

std::vector<MeshLib::Element*> GocadSGridReader::getFaceSetElements() const
{
	std::vector<MeshLib::Element*> elements;
	std::size_t face_set_id(0);
	for (auto face_set_it = _face_sets.begin(); face_set_it != _face_sets.end(); face_set_it++) {
		for (auto face_it = face_set_it->_face_pos_and_dir.begin();
				face_it != face_set_it->_face_pos_and_dir.end(); face_it++) {
			const std::size_t i((*face_it)[0]);
			const std::size_t j((*face_it)[1]);
			const std::size_t k((*face_it)[2]);

			std::array<MeshLib::Node*, 8> element_nodes;

			element_nodes[0] = _nodes[_index_calculator(i, j, k)];
			element_nodes[1] = _nodes[_index_calculator(i + 1, j, k)];
			element_nodes[2] = _nodes[_index_calculator(i + 1, j + 1, k)];
			element_nodes[3] = _nodes[_index_calculator(i, j + 1, k)];
			element_nodes[4] = _nodes[_index_calculator(i, j, k + 1)];
			element_nodes[5] = _nodes[_index_calculator(i + 1, j, k + 1)];
			element_nodes[6] = _nodes[_index_calculator(i + 1, j + 1, k + 1)];
			element_nodes[7] = _nodes[_index_calculator(i, j + 1, k + 1)];
			elements.push_back(new MeshLib::Hex(element_nodes, face_set_id));
		}
		face_set_id++;
	}

	return elements;
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

template <typename T>
T swapEndianness(T const& v)
{
	union
	{
		T v;
		char c[4];
	} a, b;

	a.v = v;
	for (unsigned short i = 0; i < sizeof(T); i++)
		b.c[i] = a.c[sizeof(T) - i - 1];

	return b.v;
}

template <typename T>
T readValue(std::ifstream& in)
{
	T v;
	in.read(reinterpret_cast<char*>(&v), sizeof(T));
	return swapEndianness(v);
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
		coords[k % 3] = readValue<float>(in);
		if ((k + 1) % 3 == 0)
			_nodes[k/3] = new MeshLib::Node(coords, k/3);
		k++;
	}
	if (k != n * 3 && !in.eof())
		ERR("Read different number of points. Expected %d floats, got %d.\n", n * 3, k);
}

void GocadSGridReader::mapRegionFlagsToCellProperties(std::vector<Bitset> const& rf)
{
	std::size_t const n = _index_calculator._n_cells;
	_properties.resize(n);
	std::fill(_properties.begin(), _properties.end(), -1);
	// region flags are stored in each node ijk and give the region index for the
	// ijk-th cell.
	for (std::size_t k(0); k < _index_calculator._z_dim-1; k++) {
		for (std::size_t j(0); j < _index_calculator._y_dim-1; j++) {
			for (std::size_t i(0); i < _index_calculator._x_dim-1; i++) {
				// Find layers containing regions given by bits.
				// Run over bits and push back set bits
				std::set<std::size_t> layers_set;
				for (auto r = regions.begin(); r != regions.end(); ++r)
				{
					if (rf[_index_calculator(i,j,k)].test(r->bit))
					{
						// Bit is set, find a layer.
						for (std::size_t l_id = 0; l_id < layers.size(); ++l_id)
							if (layers[l_id].hasRegion(*r))
								layers_set.insert(l_id);
					}
				}
				if (layers_set.size() != 1)
					ERR("A cell %d %d %d belongs to multiple (%d) layers.", i, j, k, layers_set.size());

				_properties[_index_calculator.getCellIdx(i,j,k)] = *layers_set.begin();
			}
		}
	}
}

template <typename T>
std::vector<T> readBinaryArray(std::string const& filename, std::size_t const n)
{
	std::ifstream in(filename.c_str());
	if (!in) {
		ERR("readBinaryArray(): Error while reading from file \"%s\".\n", filename.c_str());
		ERR("Could not open file \"%s\" for input.\n", filename.c_str());
		in.close();
		return std::vector<T>();
	}

	std::vector<T> result;
	result.reserve(n);

	while (in)
		result.push_back(readValue<T>(in));

	if (!in.eof())
	{
		ERR("readBinaryArray(): Error while reading from file \"%s\".\n", filename.c_str());
		ERR("Input stream invalid.\n");
		return std::vector<T>();
	}

	if (result.size() - 1 != n)
	{
		ERR("readBinaryArray(): Error while reading from file \"%s\".\n", filename.c_str());
		ERR("Read different number of values. Expected %d, got %d.\n", n, result.size());
	}

	return result;
}

void GocadSGridReader::readElementPropertiesBinary()
{
	_properties = readBinaryArray<float>(_properties_fname, _index_calculator._n_cells);
	if (_properties.empty())
		ERR("Reading of element properties file \"%s\" failed.", _properties_fname.c_str());
}

std::vector<int> GocadSGridReader::readFlagsBinary() const
{
	std::vector<int> result = readBinaryArray<int32_t>(
		_flags_fname,
		_index_calculator._n_nodes);

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
	if (_properties.empty()) {
		_properties.resize(_index_calculator._n_cells);
		std::fill(_properties.begin(), _properties.end(), 0);
	}
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
				_elements[cnt] = new MeshLib::Hex(element_nodes, static_cast<unsigned>(_properties[_index_calculator.getCellIdx(i,j,k)]));
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
						return (elem->getContent() < std::numeric_limits<double>::epsilon());
					}
			)
	);
	_elements.erase(new_end, _elements.end());
}

} // end namespace FileIO
