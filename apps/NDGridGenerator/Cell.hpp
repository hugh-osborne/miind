#include "Simplex.hpp"
#include <cmath>

class Cell {
public:
    std::vector<unsigned int> grid_coords;
    unsigned int num_dimensions;
    Triangulator& triangulator;
    std::vector<Simplex> simplices;
    std::map<unsigned int, std::vector<double>> hyps;

    Cell(std::vector<unsigned int> _coords, unsigned int _num_dims, std::vector<Point>& _points, Triangulator& _triangulator):
    grid_coords(_coords),
    num_dimensions(_num_dims),
    triangulator(_triangulator) {
        simplices = generateSimplices(_points);
        hyps = calculateAAHyperplanes(_points);
    }

    std::vector<Simplex> generateSimplices(std::vector<Point>& _points) {
        return triangulator.generateCellSimplices(num_dimensions, _points);
    }

    double getVolume(){
        double vol = 0.0;
        for (Simplex s : simplices)
            vol += s.getVolume();
        return vol;
    }

    std::map<unsigned int, std::vector<double>> calculateAAHyperplanes(std::vector<Point>& _points) {
        std::map<unsigned int, std::vector<double>> out;
        for(unsigned int d=0; d<num_dimensions; d++) {
            double max = _points[0].coords[d];
            double min = _points[0].coords[d];
            for(unsigned int i=1; i<_points.size(); i++) {
                if (max < _points[i].coords[d])
                    max = _points[i].coords[d];
                if (min > _points[i].coords[d])
                    min = _points[i].coords[d];
            }
            std::vector<double> pair(2);
            pair[0] = min;
            pair[1] = max;

            out[d] = pair;
        }
        return out;
    }

    double intersectsWith(Cell& other) {
        double vol_eps = 0.00000000001;
        double orig_vol = getVolume();

        std::vector<Simplex> test_simplices = simplices;
        
        for(auto const& kv : other.hyps){
            std::vector<Simplex> new_simplices_1;
            for(Simplex s : test_simplices) {
                if (s.getVolume() < vol_eps)
                    continue;

                std::vector<Simplex> st = s.intersectWithHyperplane(kv.first, kv.second[0])[0];

                for(Simplex ns : st)
                    new_simplices_1.push_back(ns);
            }
            std::vector<Simplex> new_simplices_2;
            for(Simplex s : new_simplices_1) {
                if (s.getVolume() < vol_eps)
                    continue;

                std::vector<Simplex> st = s.intersectWithHyperplane(kv.first, kv.second[1])[1];

                for(Simplex ns : st)
                    new_simplices_2.push_back(ns);
            }
            test_simplices = new_simplices_2;
        }

        double vol_prop = 0.0;
        for(Simplex s : test_simplices){
            vol_prop += s.getVolume();
        }
        return vol_prop / orig_vol;
    }
};