#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <fstream>
#include <iostream>
#include <cfloat>
#include <vector>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;


ArrayXd find_max(VectorXd mMat)
{
	VectorXd::Index maxRow;
	double max = mMat.maxCoeff(&maxRow);
	Matrix<bool, Dynamic, 1> MM = (mMat.array() == max).matrix();
	ArrayXd Index_max = ArrayXd::Zero(mMat.size());
	for (size_t i = 0; i < MM.rows(); i++)
	{
		if (MM(i))
		{
			Index_max[i] = 1.0;
		}

	}
	return Index_max;
}

MatrixXd select_M_PN(MatrixXd C, VectorXd P)
{
	MatrixXd res(int(P.sum()), int(P.sum())), res_ing(C.rows(), int(P.sum()));
	Eigen::Index j = 0;
	for (Eigen::Index i = 0; i < C.cols(); ++i)
	{
		if (P(i) > 0) res_ing.col(j++) = C.col(i);
	}
	j = 0;
	for (Eigen::Index i = 0; i < C.cols(); ++i)
	{
		if (P(i) > 0)
			res.row(j++) = res_ing.row(i);
	}

	return res;

}

MatrixXd select_V_PN(VectorXd V, VectorXd P)
{
	VectorXd res(int(P.sum()));
	Eigen::Index j = 0;
	for (Eigen::Index i = 0; i < V.size(); ++i)
	{
		if (P(i) > 0) res(j++) = V(i);
	}
	return res;

}

VectorXd lsqnonneg(MatrixXd C, VectorXd d, double tol)
/*
min_{x>0}||Cx-d||_2^2
Reference: Lawson and Hanson,
"Solving Least Squares Problems",
Prentice-Hall, 1974.
*/
{
	int n = C.rows();
	// Initialize vector of n zeros and Infs(to be used later)
	VectorXd nZeros = VectorXd::Zero(n);
	VectorXd wz = nZeros;
	//Initialize set of non - active columns to null
	VectorXd P = nZeros;
	// Initialize set of active columns to all
	// and the initial point to zeros
	VectorXd Z = VectorXd::Ones(n);
	VectorXd x = nZeros;
	VectorXd w = d - C * x;
	//Set up iteration criterion
	int outeriter = 0;
	int iter = 0;
	int itmax = 3 * n;
	// Outer loop to put variables into set to hold positive coefficients
	while (Z.sum() > 0 && (w.array() > tol).any())
	{
		outeriter += 1;
		VectorXd z = nZeros;
		wz = (wz.array()*(1-P.array())).matrix() - P* DBL_MAX;
		wz = (wz.array() *(1-Z.array()) + w.array()*Z.array()).matrix();
		ArrayXd t = find_max(wz);
		P = (P.array() * (1 - t) + t).matrix();
		Z = (Z.array()*(1 - t)).matrix();
		if (P.sum() == 1)
		{
			VectorXd::Index maxRow;
			P.maxCoeff(&maxRow);
			z(maxRow) = d(maxRow) / C(maxRow, maxRow);
		}
		if (P.sum() >= 2)
		{
			MatrixXd C_ing = select_M_PN(C, P);
			VectorXd d_ing = select_V_PN(d, P);
			VectorXd z_ing = C_ing.ldlt().solve(d_ing);
			Eigen::Index j = 0;
			for (size_t i = 0; i < P.size(); i++)
			{
				if (P(i) > 0) z(i) = z_ing(j++);
			}
		}
		while ((z.array()*P.array() + 1 - P.array() <= 0).any())
		{
			iter = iter + 1;
			if (iter > itmax)
			{
				cout << "optimfun:lsqnonneg: IterationCountExceeded" << endl;
				x = z;
				return x;
			}

			// Find indices where intermediate solution z is approximately negative
			vector<double> alpha_v;
			for (Eigen::Index i = 0; i < P.size(); i++)
			{
				if(z(i)<=0 && P(i)>0) alpha_v.push_back(x(i)/(x(i)-z(i)));
			}
			x = x + (*min_element(alpha_v.begin(), alpha_v.end())) * (z - x);
			for (size_t i = 0; i < Z.size(); i++)
			{
				if (abs(x(i)) < tol && P(i)>0)
				{
					Z(i) = 1;
				}
			}

			P = (1 - Z.array()).matrix();
			z = nZeros;
			if (P.sum() == 1)
			{
				VectorXd::Index maxRow;
				P.maxCoeff(&maxRow);
				z(maxRow) = d(maxRow) / C(maxRow, maxRow);
			}
			if (P.sum() >= 2)
			{
				MatrixXd C_ing = select_M_PN(C, P);
				VectorXd d_ing = select_V_PN(d, P);
				VectorXd z_ing = C_ing.ldlt().solve(d_ing);
				Eigen::Index j = 0;
				for (size_t i = 0; i < P.size(); i++)
				{
					if (P(i) > 0) z(i) = z_ing(j++);
				}
			}
		}
		x = z;
		w = d - C*x;
	}
	return x;
}

MatrixXd SNMF_inner(MatrixXd C_mat, MatrixXd B_mat, double tol)
{
	int n = B_mat.cols();
	for (Eigen::Index i = 0; i < n; i++)
	{
		B_mat.col(i) = lsqnonneg(C_mat, B_mat.col(i), tol);
	}
	return B_mat;
}

namespace py = pybind11;
PYBIND11_MODULE(libSNMF_inner,m)
{
m.doc() = "SNMF_inner iteration";
m.def("SNMF_inner", &SNMF_inner);
}