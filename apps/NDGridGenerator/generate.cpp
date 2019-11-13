#include "Grid.hpp"

int main() {
	std::vector<double> base = {-15.0,-15.0,-2.0};
	std::vector<double> dims = {20.0,20.0,8.0};
	std::vector<unsigned int> res = {100,100,100};
	double threshold = 1.99;
	double reset_v = -1.99;
	Grid g(base, dims, res, threshold, reset_v, 0.1);

	g.generateModelFile("hindmarsh_rose");
	g.generateTMatFileBatched("hindmarsh_rose");
	return 0;
}
