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

namespace MathLib
{
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

void applyKnownSolution(LisMatrix &A, LisVector &b,
                        const std::vector<std::size_t> &vec_knownX_id,
                        const std::vector<double> &vec_knownX_x,
                        double penalty_scaling)
{
	const std::size_t n_bc = vec_knownX_id.size();

	for (std::size_t i_bc = 0; i_bc < n_bc; i_bc++)
		applyKnownSolution(A, b, vec_knownX_id[i_bc], vec_knownX_x[i_bc]);

}

}  // MathLib
