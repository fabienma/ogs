/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ElementCoordinatesMappingLocal.h"

#include <limits>
#include <cassert>

#include "MathLib/MathTools.h"
#include "MeshLib/Node.h"

namespace MeshLib
{

ElementCoordinatesMappingLocal::ElementCoordinatesMappingLocal(const Element* e, const CoordinateSystem &coordinate_system)
: _pt_translate(std::array<double,3>({{0,0,0}}))
{
    //assert (e->getDimension() <= coordinate_system.getDimension());

    // set initial coordinates
    for(size_t i = 0; i < e->getNNodes(); i++)
        _point_vec.push_back(MathLib::Point3d(e->getNode(i)->getCoords()));

    //flip(coordinate_system, _point_vec);
    if (e->getDimension() < coordinate_system.getDimension()) {
        translate(_point_vec, _point_vec[0]);
        rotate(*e, coordinate_system, _point_vec);
    }
}

void ElementCoordinatesMappingLocal::flip(const CoordinateSystem &coordinate_system,
                                          std::vector<MathLib::Point3d> &vec_pt)
{
    switch(coordinate_system.getType())
    {
    case CoordinateSystemType::Y:
        {
            for(size_t i = 0; i < vec_pt.size(); i++)
                std::swap(vec_pt[i][0], vec_pt[i][1]);
        }
        break;
    case CoordinateSystemType::Z:
        {
            for(size_t i = 0; i < vec_pt.size(); i++)
                std::swap(vec_pt[i][0], vec_pt[i][2]);
        }
        break;
    case CoordinateSystemType::XZ:
        {
            for(size_t i = 0; i < vec_pt.size(); i++)
                std::swap(vec_pt[i][1], vec_pt[i][2]);
        }
        break;
    case CoordinateSystemType::YZ:
        {
            for(size_t i = 0; i < vec_pt.size(); i++)
            {
                double tmp_x = vec_pt[i][0];
                vec_pt[i][0] = vec_pt[i][1];
                vec_pt[i][1] = vec_pt[i][2];
                vec_pt[i][2] = tmp_x;
            }
        }
        break;
    default:
        break;
    }
}

void ElementCoordinatesMappingLocal::translate(std::vector<MathLib::Point3d> &vec_pt, const MathLib::Point3d origin)
{
    for (std::size_t i=0; i<vec_pt.size(); ++i) {
        vec_pt[i] -= origin;
    }
    _pt_translate = origin;
}

void ElementCoordinatesMappingLocal::rotate(const Element &ele, const CoordinateSystem &coordinate_system, std::vector<MathLib::Point3d> &vec_pt)
{
    // compute a rotation matrix
    getRotationMatrixToOriginal(ele, coordinate_system, vec_pt, _matR2original);

    // rotate the point coordinates
    const std::size_t global_dim = coordinate_system.getDimension();
    double const* const coords_node_0 (vec_pt[0].getCoords());
    Eigen::VectorXd dx = Eigen::VectorXd::Zero(global_dim);
    Eigen::VectorXd x_new = Eigen::VectorXd::Zero(3);
    for(std::size_t i = 0; i < ele.getNNodes(); i++)
    {
        double const* const coords_node_i (vec_pt[i].getCoords());
        for (std::size_t j=0; j<global_dim; j++)
            dx[j] = (coords_node_i[j] - coords_node_0[j]);

        x_new.head(global_dim) = _matR2original.transpose() * dx;
        _point_vec[i] = MathLib::Point3d(x_new.data());
    }
};

// x=Rx' where x is original coordinates and x' is local coordinates
void ElementCoordinatesMappingLocal::getRotationMatrixToOriginal(
    const Element &e,
    const CoordinateSystem &coordinate_system,
    const std::vector<MathLib::Point3d> &vec_pt,
    EMatrix &matR)
{
    const std::size_t global_dim = coordinate_system.getDimension();

    matR = EMatrix::Zero(global_dim, global_dim);
    if (global_dim == e.getDimension()) {
        matR = EMatrix::Identity(global_dim, global_dim);
    } else if (global_dim == 2 && e.getDimension() == 1) {
        MathLib::Vector3 xx(vec_pt[0], vec_pt[1]);
        xx[2] = 0.0;
        xx.normalize();
        double cos_theta = xx[0];
        double sin_theta = xx[1];
        matR(0,0) = matR(1,1) = cos_theta;
        matR(0,1) = - sin_theta;
        matR(1,0) = sin_theta;
    } else if (global_dim == 3 && e.getDimension() == 2) {
        // x"_vec
        MathLib::Vector3 xx(vec_pt[0], vec_pt[1]);
        xx.normalize();
        // a vector on the plane
        MathLib::Vector3 yy(vec_pt[1], vec_pt[2]);
        // z"_vec. off plane
        MathLib::Vector3 zz(MathLib::crossProduct(xx, yy));
        zz.normalize();
        // y"_vec
        yy = MathLib::crossProduct(zz, xx);
        yy.normalize();

        for (size_t i=0; i<global_dim; ++i) {
            matR(i, 0) = xx[i];
            matR(i, 1) = yy[i];
            matR(i, 2) = zz[i];
        }
    } else if (global_dim == 3 && e.getDimension() == 1) {
        // x"_vec
        MathLib::Vector3 xx(vec_pt[0], vec_pt[1]);
        xx.normalize();
        // an arbitrary vector
        MathLib::Vector3 yy(vec_pt[0], vec_pt[1]);
        for (size_t i = 0; i < 3; i++)
            yy[i] = 0.0;

        if (fabs(xx[0]) > 0.0 && fabs(xx[1]) + fabs(xx[2]) < std::numeric_limits<double>::epsilon()) {
            yy[2] = 1.0;
        } else if (fabs(xx[1]) > 0.0 && fabs(xx[0]) + fabs(xx[2]) < std::numeric_limits<double>::epsilon()) {
            yy[0] = 1.0;
        } else if (fabs(xx[2]) > 0.0 && fabs(xx[0]) + fabs(xx[1]) < std::numeric_limits<double>::epsilon()) {
            yy[1] = 1.0;
        } else {
            for (size_t i = 0; i < 3; i++) {
                if (fabs(xx[i]) > 0.0) {
                    yy[i] = -xx[i];
                    break;
                }
            }
        }
        // z"_vec
        MathLib::Vector3 zz(MathLib::crossProduct(xx, yy));
        zz.normalize();
        // y"_vec
        yy = MathLib::crossProduct(zz, xx);
        yy.normalize();

        for (size_t i=0; i<global_dim; ++i) {
            matR(i, 0) = xx[i];
            matR(i, 1) = yy[i];
            matR(i, 2) = zz[i];
        }
    }

}

} // MeshLib
