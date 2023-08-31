#include <lr_interpolation.h>
#include <random>
#include <vector>
#include <algorithm>
#include <limits>
void main() {
	const double infinity = std::numeric_limits<double>::infinity();
	std::vector<Vec3<double>> source;
	std::vector<Vec3<double>> target;
	bool use_random = false;

	double row_min = +infinity;
	double row_max = -infinity;
	double col_min = +infinity;
	double col_max = -infinity;

	if (use_random) {
		std::default_random_engine generator;
		std::uniform_real_distribution<float> dist_pos(-100, 100);
		std::gamma_distribution<float> dist_value(3.5f, 1.0f);
		//  setup sampling function
		//auto data_point_random = [&dist_pos, &dist_value, &generator]() -> std::pair<Vec2<float>, float>
		//{
		//	float x, y, value;

		//	x = dist_pos(generator);
		//	y = dist_pos(generator);

		//	//  sample theta
		//	value = dist_value(generator);

		//	//for (int i = 0; i < 100; i++) {
		//	//	auto data_point = data_point_random();
		//	//}

		//};
	}
	else {
		std::ifstream dataFile(".\\ZoneA.csv", std::ifstream::in);

		if (dataFile.good())
		{
			std::string line;
			std::getline(dataFile, line);
			while (!dataFile.eof())
			{
				std::getline(dataFile, line);
				std::stringstream lineStream(line);
				std::string cell;
				unsigned int column = 0;
				double x = 0.0;
				double y = 0.0;
				double value = 0.0;
				while (std::getline(lineStream, cell, ',') && column < 4)
				{
					if (column == 0)
					{
						x = std::atof(cell.c_str());
					}
					else if (column == 1)
					{
						y = std::atof(cell.c_str());
					}
					else if (column == 3)
					{
						value = std::atof(cell.c_str());

						source.push_back(Vec3(x, y, value));
						row_min = std::min(x, row_min);
						row_max = std::max(x, row_max);

						col_min = std::min(y, col_min);
						col_max = std::max(y, col_max);

					}

					column++;
				}
			}
		}
	}

#define Grid_cells 10
	target.resize(Grid_cells * Grid_cells);
	for (int i = 0; i < Grid_cells; i++) {
		for (int j = 0; j < Grid_cells; j++) {
			auto x = row_min + (row_max - row_min) * i / Grid_cells;
			auto y = col_min + (col_max - col_min) * j / Grid_cells;
			target[i * Grid_cells + j].data[0] = x;
			target[i * Grid_cells + j].data[1] = y;
		}
	}


	std::list<Vec3<double>> source_array(source.begin(), source.end());
	std::list<Vec3<double>> target_array(target.begin(), target.end());
	spatial_interpolation(source_array, target_array);
}