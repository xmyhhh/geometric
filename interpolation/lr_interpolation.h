#ifndef LR_INTERPOLATION
#define LR_INTERPOLATION
#include <list>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include "kriging.h"
#include <limits>
const double infinity = std::numeric_limits<double>::infinity();

typedef enum {
	Kriging_Interpolation = 0,
	Liner_Interpolation = 1,
	IDW_Interpolation = 2,
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

	for (auto iter_t = targetDataPoint.begin(); iter_t != targetDataPoint.end(); iter_t++) {
		//Step 1: calculate d_i

		std::vector<T> vD_i;
		vD_i.resize(sourceDataPoint.size());
		for (auto iter_s = sourceDataPoint.begin(); iter_s != sourceDataPoint.end(); iter_s++) {
			int index = std::distance(sourceDataPoint.begin(), iter_s);
			vD_i[index] = std::sqrt(std::pow((*iter_s).data[0] - (*iter_t).data[0], 2) + std::pow((*iter_s).data[1] - (*iter_t).data[1], 2)) + 0.00001;
		}


		T vD_i_inverse_sum = 0.0;
		for (auto& vD_i_item : vD_i) {
			vD_i_inverse_sum += 1.0 / (vD_i_item);

			std::cout << vD_i_item << " _ " << vD_i_inverse_sum << std::endl;
		}

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
void spatial_interpolation(std::list<Vec3<T>>& sourceDataPoint, std::list<Vec3<T>>& targetDataPoint, InterpolationMethod method = IDW_Interpolation) {

	//if (method == Liner_Interpolation)
	//	return interpolation_liner(position, value, interpolation_value);
	if (method == Kriging_Interpolation)
		return spatial_interpolation_kriging(sourceDataPoint, targetDataPoint);
	if (method == IDW_Interpolation)
		return spatial_interpolation_idw(sourceDataPoint, targetDataPoint);
}
#endif