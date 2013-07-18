/**
 * \file
 * \author Karsten Rink
 * \date   2010-03-19
 * \brief  Implementation of the VtkMeshSource class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkMeshSource.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "Elements/Element.h"
#include "Mesh.h"
#include "Node.h"
#include "VtkColorLookupTable.h"

#include "Color.h"

// ** VTK INCLUDES **
#include "vtkObjectFactory.h"
#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>

// OGS Cell Types
#include <vtkHexahedron.h>
#include <vtkLine.h>
#include <vtkPyramid.h>
#include <vtkQuad.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkWedge.h> // == Prism

vtkStandardNewMacro(VtkMeshSource);
vtkCxxRevisionMacro(VtkMeshSource, "$Revision$");

VtkMeshSource::VtkMeshSource() :
		_grid(nullptr)
{
	_removable = false; // From VtkAlgorithmProperties
	this->SetNumberOfInputPorts(0);

	const GeoLib::Color* c = GeoLib::getRandomColor();
	vtkProperty* vtkProps = GetProperties();
	vtkProps->SetColor((*c)[0] / 255.0,(*c)[1] / 255.0,(*c)[2] / 255.0);
	delete c;
	vtkProps->SetEdgeVisibility(1);
}

VtkMeshSource::~VtkMeshSource()
{
}

void VtkMeshSource::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);

	if (_grid == nullptr)
		return;
	const std::vector<MeshLib::Node*> nodes = _grid->getNodes();
	const std::vector<MeshLib::Element*> elems = _grid->getElements();
	if (nodes.empty() || elems.empty() )
		return;

	os << indent << "== VtkMeshSource ==" << "\n";

	int i = 0;
	for (std::vector<MeshLib::Node*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it)
	{
		os << indent << "Point " << i << " (" << (*it)[0] << ", " << (*it)[1] << ", " << (*it)[2] << ")\n";
	}

	i = 0;
	for (std::vector<MeshLib::Element*>::const_iterator it = elems.begin(); it != elems.end(); ++it)
	{
		os << indent << "Element " << i << ": ";
		for (unsigned t = 0; t < (*it)->getNNodes(); ++t)
			os << (*it)->getNode(t)->getID() << " ";
		os << "\n";
	}
}

int VtkMeshSource::RequestData( vtkInformation* request,
                                vtkInformationVector** inputVector,
                                vtkInformationVector* outputVector )
{
	(void)request;
	(void)inputVector;

	if (_grid == nullptr)
		return 0;
	const std::vector<MeshLib::Node*> nodes = _grid->getNodes();
	const std::vector<MeshLib::Element*> elems = _grid->getElements();

	const size_t nPoints = _grid->getNNodes();
	const size_t nElems  = _grid->getNElements();
	if (nPoints == 0 || nElems == 0)
		return 0;

	vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);
	vtkSmartPointer<vtkUnstructuredGrid> output = vtkUnstructuredGrid::SafeDownCast(
	        outInfo->Get(vtkDataObject::DATA_OBJECT()));
	output->Allocate(nElems);

	if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) > 0)
		return 1;

	// Insert grid points
	vtkSmartPointer<vtkPoints> gridPoints = vtkSmartPointer<vtkPoints>::New();
	gridPoints->SetNumberOfPoints(nPoints);
	// Generate mesh nodes
	for (unsigned i = 0; i < nPoints; ++i)
		gridPoints->SetPoint(i, (*nodes[i])[0], (*nodes[i])[1], (*nodes[i])[2]);

	// Generate mesh elements
	for (unsigned i = 0; i < nElems; ++i)
	{
		int type(0);
		const MeshLib::Element* elem (elems[i]);

		vtkIdList* point_ids = vtkIdList::New();
		const unsigned nElemNodes (elem->getNNodes());
		point_ids->SetNumberOfIds(nElemNodes);
		for (unsigned j = 0; j < nElemNodes; ++j)
			point_ids->SetId(j, elem->getNode(j)->getID());

		switch (elem->getGeomType())
		{
		case MeshElemType::LINE:
			type = 3;
			break;
		case MeshElemType::TRIANGLE:
			type = 5;
			break;
		case MeshElemType::QUAD:
			type = 9;
			break;
		case MeshElemType::HEXAHEDRON:
			type = 12;
			break;
		case MeshElemType::TETRAHEDRON:
			type = 10;
			break;
		case MeshElemType::PRISM:
			type = 13;
			for (unsigned i=0; i<3; ++i)
			{
				const unsigned prism_swap_id = point_ids->GetId(i);
				point_ids->SetId(i, point_ids->GetId(i+3));
				point_ids->SetId(i+3, prism_swap_id);
			}
			break;
		case MeshElemType::PYRAMID:
			type = 14;
			break;
		default: // if none of the above can be applied
			ERR("VtkMeshSource::RequestData(): Unknown element type \"%s\".",
					MeshElemType2String(elem->getGeomType()).c_str());
			return 0;
		}

		materialIDs->InsertValue(i, elem->getValue());
		vtkIdList* point_ids = vtkIdList::New();
		const unsigned nElemNodes (elem->getNNodes());
		point_ids->SetNumberOfIds(nElemNodes);
		for (unsigned j = 0; j < nElemNodes; ++j)
			point_ids->SetId(j, elem->getNode(j)->getID());

		output->InsertNextCell(type, point_ids);
	}

	output->SetPoints(gridPoints);

	std::vector<std::string> prop_names(_grid->getPropertyVecNames(false));
	for (auto it = prop_names.begin(); it != prop_names.end(); it++) {
		boost::optional<std::vector<unsigned> const&> grid_prop(_grid->getUnsignedPropertyVec(*it));
		if (grid_prop) {
			vtkSmartPointer<vtkIntArray> prop = vtkSmartPointer<vtkIntArray>::New();
			prop->SetName(it->c_str());
			const std::size_t tuple_size((*grid_prop).size() / nElems);
			prop->SetNumberOfComponents(tuple_size);
			prop->SetNumberOfTuples(nElems);

			for (unsigned i = 0; i < nElems*tuple_size; ++i) {
				prop->InsertValue(i, (*grid_prop)[i]);
			}

			output->GetCellData()->AddArray(prop);
			output->GetCellData()->SetActiveAttribute(it->c_str(), vtkDataSetAttributes::SCALARS);
		}
	}

	prop_names = _grid->getPropertyVecNames(true);
	for (auto it = prop_names.begin(); it != prop_names.end(); it++) {
		boost::optional<std::vector<double> const&> grid_prop(_grid->getDoublePropertyVec(*it));
		if (grid_prop) {
			vtkSmartPointer<vtkDoubleArray> prop = vtkSmartPointer<vtkDoubleArray>::New();
			prop->SetName(it->c_str());
			const std::size_t tuple_size((*grid_prop).size() / nElems);
			prop->SetNumberOfComponents(tuple_size);
			prop->SetNumberOfTuples(nElems);

			for (unsigned i = 0; i < nElems * tuple_size; ++i) {
				prop->InsertValue(i, (*grid_prop)[i]);
			}
			output->GetCellData()->AddArray(prop);
		}
	}

	output->Squeeze();

	return 1;
}

void VtkMeshSource::SetUserProperty( QString name, QVariant value )
{
	VtkAlgorithmProperties::SetUserProperty(name, value);

	(*_algorithmUserProperties)[name] = value;
}
