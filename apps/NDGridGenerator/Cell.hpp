#include "Simplex.hpp"

class Cell {
public:
    std::vector<unsigned int> grid_coords;
    unsigned int num_dimesions;
    Triangulator& triangulator;
    std::vector<Point*> points;
    std::vector<Simplex> simplices;
    std::map<unsigned int, std::vector<double>> hyps;

    Cell(std::vector<unsigned int> _coords, unsigned int _num_dims, std::vector<Point*> _points, Triangulator& _triangulator):
    grid_coords(_coords),
    num_dimesions(_num_dims),
    points(_points),
    triangulator(_triangulator) {
        simplices = generateSimplices();
        hyps = calculateAAHyperplanes();
    }

    std::vector<Simplex> generateSimplices() {
        // Based on js code by Mikola Lysenko (2014)
        // https://github.com/mikolalysenko/triangulate-hypercube

        unsigned int dfac = 1;
        for (unsigned int d=0; d<num_dimesions;d++)
            dfac *= (d+1);

        simplices = std::vector<Simplex>(dfac);
        for(unsigned int i=0; i<dfac; i++) {
            
        }

    }
};