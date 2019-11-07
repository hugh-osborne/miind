#include "Grid.hpp"

int main() {
	std::vector<double> base = {-80.0,-0.2,-0.2};
	std::vector<double> dims = {150.0,1.4,1.4};
	std::vector<unsigned int> res = {5,5,5};
	double threshold = 29.0;
	Grid g(base, dims, res, threshold, 0.1);
	g.calculateTransitionMatrix();
	return 0;
}
