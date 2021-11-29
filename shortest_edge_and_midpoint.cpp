// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "shortest_edge_and_midpoint.h"
#include <Eigen/LU>
#include <iostream>
IGL_INLINE void igl::shortest_edge_and_midpoint(
	const int e,
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& /*F*/,
	const Eigen::MatrixXi& E,
	const Eigen::VectorXi& /*EMAP*/,
	const Eigen::MatrixXi& /*EF*/,
	const Eigen::MatrixXi& /*EI*/,
	double& cost,
	Eigen::RowVectorXd& p)
{
	cost = (V.row(E(e, 0)) - V.row(E(e, 1))).norm();
	p = 0.5 * (V.row(E(e, 0)) + V.row(E(e, 1)));
}


IGL_INLINE void lowest_error_cost(
	const int e,
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& /*F*/,
	const Eigen::MatrixXi& E,
	const Eigen::VectorXi& /*EMAP*/,
	const Eigen::MatrixXi& /*EF*/,
	const Eigen::MatrixXi& /*EI*/,
	double& cost,
	Eigen::RowVectorXd& p,
	std::vector<Eigen::Matrix4d> Qs)
{
	using namespace std;
	Eigen::Matrix4d q1andq2 = (Qs.at(E(e, 0)) + Qs.at(E(e, 1)));
	q1andq2.row(3) = Eigen::RowVector4d(0, 0, 0, 1);
	if (q1andq2.determinant() != 0 && !isnan(q1andq2.determinant())) {
		p = q1andq2.inverse() * Eigen::Vector4d(0, 0, 0, 1);

	}
	else {
		p = 0.5 * (V.row(E(e, 0)) + V.row(E(e, 1)));
		Eigen::Vector4d temp = { p[0], p[1], p[2], 1 };
		p = temp;
	}
	cost = (p * (Qs.at(E(e, 0)) + Qs.at(E(e, 1))) * p.transpose());
}
