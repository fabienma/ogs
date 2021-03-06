/**
* \file
* \author Lars Bilke
* \date   2014-02-26
* \brief  Unit tests for In-Situ data arrays
*
* \copyright
* Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/
#include <numeric>

#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>

#include "gtest/gtest.h"

#include "InSituLib/VtkMappedPropertyVectorTemplate.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"

// Creates a PropertyVector<double> and maps it into a vtkDataArray-equivalent
TEST(InSituLibMappedPropertyVector, Double)
{
	const size_t mesh_size = 5;
	const double length = 1.0;

	MeshLib::Mesh* mesh = MeshLib::MeshGenerator::generateRegularHexMesh(length, mesh_size);

	ASSERT_TRUE(mesh != nullptr);
	const std::size_t number_of_tuples(mesh_size*mesh_size*mesh_size);

	std::string const prop_name("TestProperty");
	boost::optional<MeshLib::PropertyVector<double> &> double_properties(
		mesh->getProperties().createNewPropertyVector<double>(prop_name,
			MeshLib::MeshItemType::Cell));
	(*double_properties).resize(number_of_tuples);
	std::iota((*double_properties).begin(), (*double_properties).end(), 1);

	vtkNew<InSituLib::VtkMappedPropertyVectorTemplate<double> > dataArray;
	dataArray->SetPropertyVector(*double_properties);

	ASSERT_EQ(dataArray->GetNumberOfComponents(), 1);
	ASSERT_EQ(dataArray->GetNumberOfTuples(), number_of_tuples);

	ASSERT_EQ(dataArray->GetValueReference(0), 1.0);
	double* range = dataArray->GetRange(0);
	ASSERT_EQ(range[0], 1.0);
	ASSERT_EQ(range[1], 1.0 + mesh->getNElements() - 1.0);

	delete mesh;
}

// Creates a PropertyVector<int> and maps it into a vtkDataArray-equivalent
TEST(InSituLibMappedPropertyVector, Int)
{
	const size_t mesh_size = 5;
	const double length = 1.0;

	MeshLib::Mesh* mesh = MeshLib::MeshGenerator::generateRegularHexMesh(length, mesh_size);

	ASSERT_TRUE(mesh != nullptr);
	const std::size_t number_of_tuples(mesh_size*mesh_size*mesh_size);

	std::string const prop_name("TestProperty");
	boost::optional<MeshLib::PropertyVector<int> &> properties(
		mesh->getProperties().createNewPropertyVector<int>(prop_name,
			MeshLib::MeshItemType::Cell));
	(*properties).resize(number_of_tuples);
	std::iota((*properties).begin(), (*properties).end(), -5);

	vtkNew<InSituLib::VtkMappedPropertyVectorTemplate<int> > dataArray;
	dataArray->SetPropertyVector(*properties);

	ASSERT_EQ(dataArray->GetNumberOfComponents(), 1);
	ASSERT_EQ(dataArray->GetNumberOfTuples(), number_of_tuples);

	ASSERT_EQ(dataArray->GetValueReference(0), -5);
	double* range = dataArray->GetRange(0);
	ASSERT_EQ(range[0], -5);
	ASSERT_EQ(range[1], -5 + mesh->getNElements() - 1);

	delete mesh;
}

// Creates a PropertyVector<unsigned> and maps it into a vtkDataArray-equivalent
TEST(InSituLibMappedPropertyVector, Unsigned)
{
	const size_t mesh_size = 5;
	const double length = 1.0;

	MeshLib::Mesh* mesh = MeshLib::MeshGenerator::generateRegularHexMesh(length, mesh_size);

	ASSERT_TRUE(mesh != nullptr);
	const std::size_t number_of_tuples(mesh_size*mesh_size*mesh_size);

	std::string const prop_name("TestProperty");
	boost::optional<MeshLib::PropertyVector<unsigned> &> properties(
		mesh->getProperties().createNewPropertyVector<unsigned>(prop_name,
														   MeshLib::MeshItemType::Cell));
	(*properties).resize(number_of_tuples);
	std::iota((*properties).begin(), (*properties).end(), 0);

	vtkNew<InSituLib::VtkMappedPropertyVectorTemplate<unsigned> > dataArray;
	dataArray->SetPropertyVector(*properties);

	ASSERT_EQ(dataArray->GetNumberOfComponents(), 1);
	ASSERT_EQ(dataArray->GetNumberOfTuples(), number_of_tuples);

	ASSERT_EQ(dataArray->GetValueReference(0), 0);
	double* range = dataArray->GetRange(0);
	ASSERT_EQ(range[0], 0);
	ASSERT_EQ(range[1], 0 + mesh->getNElements() - 1);

	delete mesh;
}
