#include "Grid.hpp"

int main() {
	std::vector<double> base = {2.5,-13.0,-3.0};
	std::vector<double> dims = {1.0,16.0,6.0};
	std::vector<unsigned int> res = {100,100,100};
	double threshold = 1.99;
	double reset_v = -1.99;
	Grid g(base, dims, res, threshold, reset_v, 0.1);

	g.generateModelFile("hindmarsh_rose");
	g.generateTMatFileBatched("hindmarsh_rose");
	return 0;
}
