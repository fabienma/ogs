/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation of Lis utility functions.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LisTools.h"

#include <cassert>

#include "logog/include/logog.hpp"

#include "LisMatrix.h"
#include "LisVector.h"

#include "BaseLib/quicksort.h"
#include "MathLib/LinAlg/Sparse/CRSMatrix.h"

namespace MathLib
{

namespace detail
{
MathLib::CRSMatrix<double, unsigned>* lis2crs(LIS_MATRIX &A)
{
	unsigned const n(A->n);
	unsigned *iA(new unsigned[n+1]);

	iA[0] = 0;
	for (LIS_INT k=1; k<n+1; ++k) {
		iA[k] = iA[k-1] + A->w_row[k-1 - A->is];
	}

	unsigned *jA(new unsigned[iA[n]]);
	double *entries(new double[iA[n]]);
	for (unsigned r(0); r<n; ++r) {
		unsigned const beg_idx(iA[r]);
		unsigned const end_idx(iA[r+1]);
		for (unsigned j(beg_idx); j<end_idx; ++j) {
			jA[j] = A->w_index[r-A->is][j-beg_idx];
			entries[j] = A->w_value[r-A->is][j-beg_idx];
		}
	}

	for (unsigned r(0); r<n; ++r) {
		unsigned const beg_idx(iA[r]);
		unsigned const end_idx(iA[r+1]);
		// sort the column entries of the row
		BaseLib::quicksort(jA, beg_idx, end_idx, entries);
	}

	return new MathLib::CRSMatrix<double,unsigned>(A->n, iA, jA, entries);
}
} // end namespace detail

void applyKnownSolution(LisMatrix &eqsA, LisVector &eqsRHS,
	const std::vector<std::size_t> &input_rows,
	const std::vector<double> &input_vals)
{
	LIS_MATRIX &A = eqsA.getRawMatrix();

	// unfortunatly the input is not sorted => copy and sort
	std::vector<std::size_t> rows(input_rows);
	std::vector<double> vals(input_vals);
	BaseLib::quicksort(rows, 0, rows.size(), vals);

	MathLib::CRSMatrix<double,unsigned> *crs_mat(MathLib::detail::lis2crs(A));
	MathLib::CRSMatrix<double,unsigned> *crs_mat_t(crs_mat->getTranspose());

	unsigned const*const iAt(crs_mat_t->getRowPtrArray());
	unsigned const*const jAt(crs_mat_t->getColIdxArray());
	double * At(const_cast<double *>(crs_mat_t->getEntryArray()));

	// b_i -= A(i,k)*val, i!=k => b_i -= A(k,i)^T * val
	// set A^T(k,i) = 0, i!=k (i.e. set column entries of original matrix to
	// zero)
	for (std::size_t r(0); r<rows.size(); ++r) {
		auto const row = rows[r];
		auto const val = vals[r];
		for (unsigned j(iAt[row]); j<iAt[row+1]; ++j) {
			if (jAt[j] == row) // skip diagonal entry
				continue;
			eqsRHS.add(jAt[j], -At[j] * val);
			At[j] = 0.0;
		}
	}

	delete crs_mat;
	crs_mat = crs_mat_t->getTranspose();
	unsigned const*const iA(crs_mat->getRowPtrArray());
	unsigned const*const jA(crs_mat->getColIdxArray());
	double * entries(const_cast<double*>(crs_mat->getEntryArray()));

	// set row entries, except the diagonal entry, to zero
	for (std::size_t r(0); r<rows.size(); ++r) {
		auto const row = rows[r];
		for (unsigned j(iA[row]); j<iA[row+1]; ++j) {
			if (jA[j] == row)
				entries[j] = 1.0;
			else
				entries[j] = 0.0;
		}
	}

	// reset the entries in the lis matrix
	unsigned cnt(0);
	for (LIS_INT row_i = 0; row_i < A->n; ++row_i) {
		for (LIS_INT j = 0; j < A->w_row[row_i - A->is]; ++j) {
			A->w_index[row_i-A->is][j] = jA[cnt];
			A->w_value[row_i-A->is][j] = entries[cnt];
			cnt++;
		}
	}
	delete crs_mat;

	// b_k = val
	for (std::size_t k(0); k<rows.size(); ++k) {
		eqsRHS.set(rows[k], vals[k]);
	}
}

/*
static void applyKnownSolution(LisMatrix &eqsA, LisVector &eqsRHS, LIS_INT row,
                               double val)
{
	LIS_MATRIX &A = eqsA.getRawMatrix();

	// set all entries of the row of the matrix to zero except the diagonal
	// entry that is set to one
	eqsA.setValue(row, row, 1);
	for (LIS_INT k = 0; k < A->w_row[row - A->is]; ++k)
	{
		assert(k < A->w_row[row - A->is]);
		if (A->w_index[row - A->is][k] != row)
			A->w_value[row - A->is][k] = 0;
	}

	// b_i -= A(i,k)*val, i!=k and set the entries of the k-th column of the
	// matrix to zero except the diagonal entry A(k, k)
	for (LIS_INT row_i = 0; row_i < A->n; ++row_i)
	{
		if (row_i == row)
			continue;
		for (LIS_INT k = 0; k < A->w_row[row_i - A->is]; ++k)
		{
			LIS_INT column_j = A->w_index[row_i - A->is][k];
			if (column_j == row)
			{
				double const A_i_j = A->w_value[row_i - A->is][k];
				eqsRHS.add(row_i, -A_i_j * val);
				A->w_value[row_i - A->is][k] = 0;
			}
		}
	}

	// b_k = val
	eqsRHS.set(row, val);
}
*/

/*
void applyKnownSolution(LisMatrix &A, LisVector &b,
                        const std::vector<std::size_t> &vec_knownX_id,
                        const std::vector<double> &vec_knownX_x,
                        double penalty_scaling)
{
	const std::size_t n_bc = vec_knownX_id.size();

	for (std::size_t i_bc = 0; i_bc < n_bc; i_bc++)
		applyKnownSolution(A, b, vec_knownX_id[i_bc], vec_knownX_x[i_bc]);

}
*/
}  // MathLib
