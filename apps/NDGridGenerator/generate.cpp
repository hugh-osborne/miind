#include "Grid.hpp"

int main() {
	std::vector<double> base = {-80.0,-0.2,-0.2};
	std::vector<double> dims = {150.0,1.4,1.4};
	std::vector<unsigned int> res = {100,100,100};
	double threshold = 29.0;
	double reset_v = -70;
	Grid g(base, dims, res, threshold, reset_v, 0.1);
	g.generateModelFile("test");
	// g.generateTMatFile("test");
	// g.generateTMatFileLowMemory("test");
	g.generateTMatFileBatched("test");
	return 0;
}
