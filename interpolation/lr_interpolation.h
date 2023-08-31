#ifndef LR_INTERPOLATION
#define LR_INTERPOLATION
#include <list>
#include <assert.h>
#include <iostream>

#include "kriging.h"
typedef enum {
	Kriging_Interpolation = 0,
	Liner_Interpolation = 1
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

	double estimatedZ = 0.0;

	for (auto iter = targetDataPoint.begin(); iter != targetDataPoint.end(); iter++) {
		Point point((*iter).data[0], (*iter).data[1]);
		estimatedZ = kriging.SimpleKrigeForPoint(point.x, point.y, model, nugget, sill, range, choleskyDecomposition, residuals, estimatedMean);
		(*iter).data[2] = estimatedZ;
	}

	//for (unsigned int i = 0; i < rasterContext.width; i++)
	//{
	//	for (unsigned int j = 0; j < rasterContext.height; j++)
	//	{
	//		// The kriged estimate

	//		Point point = rasterContext.XYtoPoint(i, j);

	//		estimatedZ = SimpleKrigeForPoint(point.x, point.y, model, nugget, sill, range, choleskyDecomposition, residuals, estimatedMean);

	//		dataPoint_value.push_back(DataPoint(i, j, estimatedZ));
	//	}
	//}
	return;
}

template <typename T = double>
void spatial_interpolation(std::list<Vec3<T>>& sourceDataPoint, std::list<Vec3<T>>& targetDataPoint, InterpolationMethod method = Kriging_Interpolation) {

	//if (method == Liner_Interpolation)
	//	return interpolation_liner(position, value, interpolation_value);
	if (method == Kriging_Interpolation)
		return spatial_interpolation_kriging(sourceDataPoint, targetDataPoint);
}
#endif