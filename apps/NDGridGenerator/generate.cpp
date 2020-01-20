#include "Grid.hpp"

int main() {
	std::vector<double> base = {-2.0,-2.0,-75};
	std::vector<double> dims = {48.0,48.0,40};
	std::vector<unsigned int> res = {50,50,200};
	double threshold = -40.4;
	double reset_v = -70.6;
	Grid g(base, dims, res, threshold, reset_v, 0.01);

	g.generateModelFile("conductanceInNdNoise", 0.001);
	g.generateTMatFileBatched("conductanceInNdNoise");
	return 0;
}
