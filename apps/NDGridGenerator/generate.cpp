#include "Grid.hpp"

int main() {
	std::vector<double> base = {-0.2,-0.2,-80.0};
	std::vector<double> dims = {1.2,1.2,150.0};
	std::vector<unsigned int> res = {50,50,50};
	double threshold = 29.0;
	double reset_v = -70;
	Grid g(base, dims, res, threshold, reset_v, 0.1);

	g.generateModelFile("rybak_3d");
	g.generateTMatFileBatched("rybak_3d");
	return 0;
}
