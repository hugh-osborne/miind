#include "Point.hpp"
#include "Cell.hpp"
#include "Simplex.hpp"
#include "Triangulator.hpp"


std::vector<Simplex> Triangulator::chooseTriangulation(unsigned int num_dimensions, std::vector<Point>& points, std::vector<unsigned int>& lower_inds, std::vector<unsigned int>& upper_inds, std::vector<unsigned int>& hyper_inds, std::vector<unsigned int>& all_inds) {
	if (lower_inds.size() == 0){
		std::vector<Simplex> out;
		std::vector<Point> ps(upper_inds.size());
		for (unsigned int i=0; i<upper_inds.size(); i++)
			ps[i] = points[upper_inds[i]];
		out.push_back(Simplex(num_dimensions, ps, *this));
		return out;
	}
	if (upper_inds.size() == 0){
		std::vector<Simplex> out;
		std::vector<Point> ps(lower_inds.size());
		for (unsigned int i=0; i<lower_inds.size(); i++)
			ps[i] = points[lower_inds[i]];
		out.push_back(Simplex(num_dimensions, ps, *this));
		return out;
	}
	std::vector<std::vector<unsigned int>> tris = transitions[lower_inds.size()][upper_inds.size()][hyper_inds.size()];
	std::vector<Simplex> out;
	for (unsigned int t=0; t <tris.size(); t++) {
		std::vector<Point> ps(tris[t].size());
		for (unsigned int i=0; i<tris[t].size(); i++)
			ps[i] = points[all_inds[i]];
		out.push_back(Simplex(num_dimensions, ps, *this));
	}
}

std::vector<Simplex> Triangulator::generateCellSimplices(unsigned int num_dimensions, std::vector<Point>& points) {
	switch(num_dimensions) {
	case 2: {
		std::vector<Simplex> simplices;
		std::vector<Point> ps_0(3);
		ps_0[0] = points[0];
		ps_0[1] = points[1];
		ps_0[2] = points[3];
		simplices.push_back(Simplex(num_dimensions,ps_0,*this));
		std::vector<Point> ps_1(3);
		ps_1[0] = points[3];
		ps_1[1] = points[2];
		ps_1[2] = points[0];
		simplices.push_back(Simplex(num_dimensions,ps_1,*this));
		return simplices;
	}
	case 3: {
		std::vector<Simplex> simplices;
		std::vector<Point> ps_0(4);
		ps_0[0] = points[0];
		ps_0[1] = points[1];
		ps_0[2] = points[3];
		ps_0[3] = points[7];
		simplices.push_back(Simplex(num_dimensions,ps_0,*this));
		std::vector<Point> ps_1(4);
		ps_1[0] = points[7];
		ps_1[1] = points[1];
		ps_1[2] = points[5];
		ps_1[3] = points[0];
		simplices.push_back(Simplex(num_dimensions,ps_1,*this));
		std::vector<Point> ps_2(4);
		ps_2[0] = points[7];
		ps_2[1] = points[2];
		ps_2[2] = points[3];
		ps_2[3] = points[0];
		simplices.push_back(Simplex(num_dimensions,ps_2,*this));
		std::vector<Point> ps_3(4);
		ps_3[0] = points[0];
		ps_3[1] = points[2];
		ps_3[2] = points[6];
		ps_3[3] = points[7];
		simplices.push_back(Simplex(num_dimensions,ps_3,*this));
		std::vector<Point> ps_4(4);
		ps_4[0] = points[0];
		ps_4[1] = points[4];
		ps_4[2] = points[5];
		ps_4[3] = points[7];
		simplices.push_back(Simplex(num_dimensions,ps_4,*this));
		std::vector<Point> ps_5(4);
		ps_5[0] = points[7];
		ps_5[1] = points[4];
		ps_5[2] = points[6];
		ps_5[3] = points[0];
		simplices.push_back(Simplex(num_dimensions,ps_5,*this));
		return simplices;
	}
	default: {
		return std::vector<Simplex>();
	}
	}
}
