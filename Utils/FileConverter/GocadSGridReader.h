/**
 * @file GocadSGridReader.h
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

#ifndef GOCADSGRIDREADER_H_
#define GOCADSGRIDREADER_H_

// STL
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

// MeshLib
#include "Node.h"
#include "Elements/Element.h"

#include <boost/dynamic_bitset.hpp>

namespace FileIO
{

class GocadSGridReader
{
private:
	typedef boost::dynamic_bitset<> Bitset;
public:
	/**
	 * Constructor takes as argument the Gocad .sg text file.
	 * @param fname file name
	 */
	GocadSGridReader(std::string const& fname);
	~GocadSGridReader();

	std::vector<MeshLib::Node*> const& getNodes() const { return _nodes; }
	std::vector<MeshLib::Element*> const& getElements() const { return _elements; }
	std::vector<double> const& getCellProperties() const { return _properties; }

private:
	void parseDims(std::string const& line);
	void parsePointsFileName(std::string const& line);
	void parsePropertiesFileName(std::string const& line);
	void parseFlagsFileName(std::string const& line);
	void parseRegionFlagsFileName(std::string const& line);

	/**
	 * Class for calculating the index to given 3d position within the structured grid.
	 */
	class IndexCalculator
	{
	public:
		/**
		 * Constructor initializes the dimensions.
		 * @param x_dim
		 * @param y_dim
		 * @param z_dim
		 */
		IndexCalculator(std::size_t x_dim, std::size_t y_dim, std::size_t z_dim) :
				_x_dim(x_dim), _y_dim(y_dim), _z_dim(z_dim),
				_n_nodes(x_dim * y_dim * z_dim),
				_n_cells((_x_dim - 1) * (_y_dim - 1) * (_z_dim - 1))
		{}

		IndexCalculator() :
				_x_dim(0), _y_dim(0), _z_dim(0), _n_nodes(0), _n_cells(0)
		{}

		/**
		 * Get the index for the space coordinates i,j,k.
		 * @param i
		 * @param j
		 * @param k
		 * @return index within the field of nodes
		 */
		std::size_t operator()(std::size_t i, std::size_t j, std::size_t k)
		{
			const std::size_t idx(k * _x_dim * _y_dim + j * _x_dim + i);
			if (idx >= _n_nodes) {
				return std::numeric_limits < std::size_t > ::max();
			}
			return idx;
		}

		std::size_t getCellIdx(std::size_t i, std::size_t j, std::size_t k)
		{
			const std::size_t idx(k * (_x_dim-1) * (_y_dim-1) + j * (_x_dim-1) + i);
			if (idx >= _n_cells) {
				return std::numeric_limits < std::size_t > ::max();
			}
			return idx;
		}

		std::size_t _x_dim;
		std::size_t _y_dim;
		std::size_t _z_dim;
		std::size_t _n_nodes;
		std::size_t _n_cells;
	};

public:
	struct Region
	{
		std::string name;
		unsigned bit;

		bool operator==(Region const& r) const { return bit == r.bit; }
	};

	// Each model layer own several regions.
	struct Layer
	{
		std::vector<Region> regions;

		bool hasRegion(Region const& r) const
		{
			return (std::find(regions.begin(), regions.end(), r) != regions.end());
		}

	};

private:
	void readNodesBinary();
	std::vector<int> readFlagsBinary() const;
	std::vector<Bitset> readRegionFlagsBinary() const;
	void makeNodesUnique();
	void readElementPropertiesBinary();
	void mapRegionFlagsToCellProperties(std::vector<Bitset> const& rf);
	void createElements();
	void readSplitNodesAndModifyElements();
	void modifyElement(std::size_t u, std::size_t v, std::size_t w, MeshLib::Node const* node2sub,
			MeshLib::Node * substitute_node);
	void removeNullVolumeElements();

	std::string const& _fname;
	std::string const _path;
	// data read from sg file
	IndexCalculator _index_calculator;
	std::string _pnts_fname;
	std::string _properties_fname;
	std::string _flags_fname;
	std::string _region_flags_fname;

	std::vector<Region> regions;
	std::vector<Layer> layers;

	// data read from binary points file
	std::vector<MeshLib::Node*> _nodes;
	std::vector<double> _properties;
	// calculated data
	std::vector<std::size_t> _node_id_map;
	std::vector<MeshLib::Element*> _elements;
}; // end class GocadSGridReader

} // end namespace FileIO

#endif /* GOCADSGRIDREADER_H_ */
