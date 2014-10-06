/**
 * \file
 * \author Karsten Rink
 * \date   2011-05-10
 * \brief  Implementation of the ElementTreeModel class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementTreeModel.h"
#include "OGSError.h"
#include "TreeItem.h"
#include "Mesh.h"
#include "Node.h"
#include "MeshInformation.h"
#include "Elements/Element.h"
#include "AABB.h"

#include "VtkMeshSource.h"

/**
 * Constructor.
 */
ElementTreeModel::ElementTreeModel( QObject* parent )
	: TreeModel(parent), _mesh_source(nullptr)
{
	QList<QVariant> rootData;
	delete _rootItem;
	rootData << "Name" << "Type" << "" << "";
	_rootItem = new TreeItem(rootData, nullptr);
}

ElementTreeModel::~ElementTreeModel()
{
}

void ElementTreeModel::setElement(vtkUnstructuredGridAlgorithm const*const grid, const unsigned elem_index)
{
	this->_mesh_source = grid;
	this->clearView();

	VtkMeshSource const*const source = dynamic_cast<VtkMeshSource const*const>(grid);

	if (!source)
		return;

	const MeshLib::Mesh* mesh = source->GetMesh();
	const MeshLib::Element* elem = mesh->getElement(elem_index);
    boost::optional<std::vector<unsigned> const&> mat_ids = mesh->getUnsignedPropertyVec("MaterialIDs");

	QList<QVariant> elemData;
	elemData << "Element " + QString::number(elem_index) << "" << "" << "";
	TreeItem* elemItem = new TreeItem(elemData, _rootItem);
	_rootItem->appendChild(elemItem);

	QList<QVariant> typeData;
	typeData << "Element Type: " << QString::fromStdString(MeshElemType2String(elem->getGeomType()));
	TreeItem* typeItem = new TreeItem(typeData, elemItem);
	elemItem->appendChild(typeItem);

	if (mat_ids) {
		QList<QVariant> materialData;
		materialData << "MaterialID: " << QString::number((*mat_ids)[elem_index]);
		TreeItem* matItem = new TreeItem(materialData, elemItem);
		elemItem->appendChild(matItem);
	}

	QList<QVariant> volData;
	volData << "Area/Volume: " <<
	QString::number(mesh->getElement(elem_index)->getContent());
	TreeItem* volItem = new TreeItem(volData, elemItem);
	elemItem->appendChild(volItem);

	QList<QVariant> nodeListData;
	nodeListData << "Nodes" << "" << "" << "";
	TreeItem* nodeListItem = new TreeItem(nodeListData, elemItem);
	elemItem->appendChild(nodeListItem);

	//const std::vector<MeshLib::Node*> nodes_vec = grid->getNodes();
	size_t nElemNodes = elem->getNNodes();
	for (size_t i = 0; i < nElemNodes; i++)
	{
		const MeshLib::Node* node = elem->getNode(i);
		QList<QVariant> nodeData;
		nodeData << "Node " + QString::number(node->getID()) <<
		QString::number((*node)[0]) << QString::number((*node)[1]) <<
		QString::number((*node)[2]);
		TreeItem* nodeItem = new TreeItem(nodeData, nodeListItem);
		nodeListItem->appendChild(nodeItem);
	}
	reset();
}

void ElementTreeModel::clearView()
{
	_rootItem->removeChildren(0, _rootItem->childCount());
	reset();
}

void ElementTreeModel::setMesh(MeshLib::Mesh const*const mesh)
{
	this->clearView();

	if (!mesh)
		return;
	
	QList<QVariant> mesh_name;
	mesh_name << "Name:" << QString::fromStdString(mesh->getName()) << "" << "" << "";
	TreeItem* name_item = new TreeItem(mesh_name, _rootItem);
	_rootItem->appendChild(name_item);

	QList<QVariant> nodes_number;
	nodes_number << "#Nodes: " << QString::number(mesh->getNNodes()) << "" << "";
	TreeItem* nodes_item = new TreeItem(nodes_number, _rootItem);
	_rootItem->appendChild(nodes_item);

	QList<QVariant> elements_number;
	elements_number << "#Elements: " << QString::number(mesh->getNElements()) << "" << "";
	TreeItem* elements_item = new TreeItem(elements_number, _rootItem);
	_rootItem->appendChild(elements_item);

	const std::array<QString, 7> n_element_names = {{ "Lines:", "Triangles:", "Quads:", "Tetrahedra:", "Hexahedra:", "Pyramids:", "Prisms:" }};
	const std::array<unsigned, 7>& n_element_types (MeshInformation::getNumberOfElementTypes(*mesh));
	for (std::size_t i=0; i<n_element_types.size(); ++i)
	{
		if (n_element_types[i])
		{
			QList<QVariant> elements_number;
			elements_number << n_element_names[i] << QString::number(n_element_types[i]) << "" << "";
			TreeItem* type_item = new TreeItem(elements_number, elements_item);
			elements_item->appendChild(type_item);
		}
	}

	QList<QVariant> bounding_box;
	bounding_box << "Bounding Box" << "" << "" << "";
	TreeItem* aabb_item = new TreeItem(bounding_box, _rootItem);
	_rootItem->appendChild(aabb_item);

	const GeoLib::AABB<MeshLib::Node> aabb (MeshInformation::getBoundingBox(*mesh));
	const MeshLib::Node& min = aabb.getMinPoint();
	const MeshLib::Node& max = aabb.getMaxPoint();

	QList<QVariant> min_aabb;
	min_aabb << "Min:" << QString::number(min[0], 'f') << QString::number(min[1], 'f') << QString::number(min[2], 'f');
	TreeItem* min_item = new TreeItem(min_aabb, aabb_item);
	aabb_item->appendChild(min_item);

	QList<QVariant> max_aabb;
	max_aabb << "Max:" << QString::number(max[0], 'f') << QString::number(max[1], 'f') << QString::number(max[2], 'f');
	TreeItem* max_item = new TreeItem(max_aabb, aabb_item);
	aabb_item->appendChild(max_item);

	QList<QVariant> edges;
	edges << "Edge Length: " << "[" + QString::number(mesh->getMinEdgeLength(), 'f') + "," << QString::number(mesh->getMaxEdgeLength(), 'f') + "]" << "";
	TreeItem* edge_item = new TreeItem(edges, _rootItem);
	_rootItem->appendChild(edge_item);

	std::pair<unsigned, unsigned> mat_bounds (MeshInformation::getValueBounds(*mesh));
	QList<QVariant> materials;
	materials << "MaterialIDs: " << "[" + QString::number(mat_bounds.first) + "," << QString::number(mat_bounds.second) + "]" << "";
	TreeItem* mat_item = new TreeItem(materials, _rootItem);
	_rootItem->appendChild(mat_item);

	reset();
	
}

