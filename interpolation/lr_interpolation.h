#ifndef LR_INTERPOLATION
#define LR_INTERPOLATION
#include <list>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include "kriging.h"
#include <limits>
#include "Eigen/Dense"

const double infinity = std::numeric_limits<double>::infinity();

typedef enum {
	Kriging_Interpolation = 0,
	Liner_Interpolation = 1,
	IDW_Interpolation = 2,
	RBF_Interpolation = 3,
} InterpolationMethod;

template <typename T = double>
struct Vec3 {
public:
	T data[3];
	Vec3() {}

	Vec3(T a, T b) {
		data[0] = a;
		data[1] = b;
	}
	Vec3(T a, T b, T c) {
		data[0] = a;
		data[1] = b;
		data[2] = c;
	}
};

template <typename T = double>
void spatial_interpolation_kriging(std::list<Vec3<T>>& sourceDataPoint, std::list<Vec3<T>>& targetDataPoint) {
	std::vector<DataPoint> dataPoint_copy{};
	dataPoint_copy.resize(sourceDataPoint.size());

	for (auto iter = sourceDataPoint.begin(); iter != sourceDataPoint.end(); iter++) {
		int i = std::distance(sourceDataPoint.begin(), iter);
		dataPoint_copy[i].x = (*iter).data[0];
		dataPoint_copy[i].y = (*iter).data[1];
		dataPoint_copy[i].value = (*iter).data[2];
	}

	Kriging kriging(dataPoint_copy);

	kriging.Initialize();
	auto model = Kriging::Spherical;
	double nugget = 0.0;
	double sill = kriging.GetEstimatedSill();
	double range = 4000;
	// Find distances between all points and calculate covariograms

	std::vector<std::vector<double>> distanceCovariogramMatrix = kriging.CalculateCovariogramMatrix(kriging.dataPoint, model, nugget, sill, range, false);

	// Decompose the covariogram matrix 

	CholeskyDecomposition choleskyDecomposition(distanceCovariogramMatrix);

	choleskyDecomposition.Decompose();

	// Estimate mean and calculate residuals

	double estimatedMean = std::accumulate(kriging.dataPoint.begin(), kriging.dataPoint.end(), 0.0, [](double sum, const DataPoint& dataPoint) { return sum + dataPoint.value; }) / kriging.dataPoint.size();

	std::vector<double> residuals(kriging.dataPoint.size());

	std::transform(kriging.dataPoint.begin(), kriging.dataPoint.end(), residuals.begin(), [&](const DataPoint& current) { return current.value - estimatedMean; });

	T estimatedZ = 0.0;

	for (auto iter = targetDataPoint.begin(); iter != targetDataPoint.end(); iter++) {
		Point point((*iter).data[0], (*iter).data[1]);
		estimatedZ = kriging.SimpleKrigeForPoint(point.x, point.y, model, nugget, sill, range, choleskyDecomposition, residuals, estimatedMean);
		(*iter).data[2] = estimatedZ;
	}

	return;
}


template <typename T = double>
void spatial_interpolation_idw(std::list<Vec3<T>>& sourceDataPoint, std::list<Vec3<T>>& targetDataPoint) {
	//该idw插值没有考虑搜索领域实现，复杂度为O（m*n）
	for (auto iter_t = targetDataPoint.begin(); iter_t != targetDataPoint.end(); iter_t++) {
		//Step 1: calculate D_i

		std::vector<T> vD_i;
		vD_i.resize(sourceDataPoint.size());
		for (auto iter_s = sourceDataPoint.begin(); iter_s != sourceDataPoint.end(); iter_s++) {
			int index = std::distance(sourceDataPoint.begin(), iter_s);
			vD_i[index] = std::sqrt(std::pow((*iter_s).data[0] - (*iter_t).data[0], 2) + std::pow((*iter_s).data[1] - (*iter_t).data[1], 2)) + 0.00001;
		}


		T vD_i_inverse_sum = 0.0;
		for (auto& vD_i_item : vD_i) {
			vD_i_inverse_sum += 1.0 / (vD_i_item);

		}
		//Step 2: calculate W_i
		std::vector<T> vW_i;
		vW_i.resize(sourceDataPoint.size());
		for (int i = 0; i < sourceDataPoint.size(); i++) {
			auto vD_i_item = vD_i[i];
			vW_i[i] = (1.0 / (vD_i_item + 0.00001)) / vD_i_inverse_sum;
		}

		T estimatedZ = 0.0;
		for (auto iter_s = sourceDataPoint.begin(); iter_s != sourceDataPoint.end(); iter_s++) {
			int index = std::distance(sourceDataPoint.begin(), iter_s);
			estimatedZ += vW_i[index] * (*iter_s).data[2];

		}
		(*iter_t).data[2] = estimatedZ;

	}
	return;



}

template <typename T = double>
void spatial_interpolation_rbf(std::list<Vec3<T>>& sourceDataPoint, std::list<Vec3<T>>& targetDataPoint) {
	//ref  https://en.wikipedia.org/wiki/Radial_basis_function_interpolation

	auto inverse_multiquadric = [](T r) {
		double epsilon = 0.198992436915951;
		return 1.0 / std::sqrt(std::pow(epsilon * r, 2.0) + 1.0);
	};

	auto multiquadric = [](T r) {
		double epsilon = 0.198992436915951;
		double tmp = (epsilon * r);
		return std::sqrt(1.0 + (tmp * tmp));
	};

	int sourceDataPoint_size = sourceDataPoint.size();
	int targetDataPoint_size = targetDataPoint.size();

	assert(sourceDataPoint_size > 0);
	assert(targetDataPoint_size > 0);

	//Step 1: calculate mat_A
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Matrix_A(sourceDataPoint_size, sourceDataPoint_size);
	for (auto iter_i = sourceDataPoint.begin(); iter_i != sourceDataPoint.end(); iter_i++) {
		int i = std::distance(sourceDataPoint.begin(), iter_i);
		for (auto iter_j = iter_i; iter_j != sourceDataPoint.end(); iter_j++) {
			int j = std::distance(sourceDataPoint.begin(), iter_j);
			T distance;
			if (i == j)
				distance = 0.0;
			else {
				distance = std::sqrt(
					std::pow((*iter_i).data[0] - (*iter_j).data[0], 2)
					+
					std::pow((*iter_i).data[1] - (*iter_j).data[1], 2)
				);

				distance = multiquadric(distance);
			}


			Matrix_A(i, j) = distance;
			Matrix_A(j, i) = distance;
		}
	}
	//Step 2: calculate mat_b
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Matrix_b(sourceDataPoint_size, 1);
	for (auto iter_i = sourceDataPoint.begin(); iter_i != sourceDataPoint.end(); iter_i++) {
		int i = std::distance(sourceDataPoint.begin(), iter_i);
		Matrix_b(i, 0) = (*iter_i).data[2];
	}


	for (int i = 0; i < Matrix_A.rows(); i++) {
		for (int j = 0; j < Matrix_A.cols(); j++) {
			std::cout << Matrix_A(i, j) << "       ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << std::endl;

	for (int i = 0; i < Matrix_b.rows(); i++) {
		for (int j = 0; j < Matrix_b.cols(); j++) {
			std::cout << Matrix_b(i, j) << "       ";
		}
		std::cout << std::endl;
	}
	//Step 2: calculate mat_x
	Eigen::VectorXd  x = Matrix_A.householderQr().solve(Matrix_b);
	std::cout << "The solution is:\n" << x << std::endl;


}

template <typename T = double>
void spatial_interpolation(std::list<Vec3<T>>& sourceDataPoint, std::list<Vec3<T>>& targetDataPoint, InterpolationMethod method = RBF_Interpolation) {

	//if (method == Liner_Interpolation)
	//	return interpolation_liner(position, value, interpolation_value);
	if (method == Kriging_Interpolation)
		return spatial_interpolation_kriging(sourceDataPoint, targetDataPoint);
	if (method == IDW_Interpolation)
		return spatial_interpolation_idw(sourceDataPoint, targetDataPoint);
	if (method == RBF_Interpolation)
		return spatial_interpolation_rbf(sourceDataPoint, targetDataPoint);
}
#endif