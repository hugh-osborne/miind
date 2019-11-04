#include "Simplex.hpp"
#include <cmath>

class Cell {
public:
    std::vector<unsigned int> grid_coords;
    unsigned int num_dimensions;
    Triangulator& triangulator;
    std::vector<Point> points;
    std::vector<Simplex> simplices;
    std::map<unsigned int, std::vector<double>> hyps;

    Cell(std::vector<unsigned int> _coords, unsigned int _num_dims, std::vector<Point>& _points, Triangulator& _triangulator):
    grid_coords(_coords),
    num_dimensions(_num_dims),
    points(_points),
    triangulator(_triangulator) {
        simplices = generateSimplices();
        hyps = calculateAAHyperplanes();
    }

    std::vector<Simplex> generateSimplices() {
        return triangulator.generateCellSimplices(num_dimensions, points);
    }

    double getVolume(){
        double vol = 0.0;
        for (Simplex s : simplices)
            vol += s.getVolume();
        return vol;
    }

    std::map<unsigned int, std::vector<double>> calculateAAHyperplanes() {
        std::map<unsigned int, std::vector<double>> out;
        for(unsigned int d=0; d<num_dimensions; d++) {
            double max = points[0].coords[d];
            double min = points[0].coords[d];
            for(unsigned int i=1; i<points.size(); i++) {
                if (max < points[i].coords[d])
                    max = points[i].coords[d];
                if (min > points[i].coords[d])
                    min = points[i].coords[d];
            }
            std::vector<double> pair(2);
            pair[0] = min;
            pair[1] = max;

            out[d] = pair;
        }
        return out;
    }

    double intersectsWith(Cell& other) {
        double vol_eps = 0.0000000000001;
        double orig_vol = getVolume();
        
        for(auto const& kv : hyps){
            std::vector<Simplex> new_simplices_1;
            for(Simplex s : simplices) {
                if (s.getVolume() < vol_eps)
                    continue;
                for(Simplex ns : s.intersectWithHyperplane(kv.first, kv.second[0])[1])
                    new_simplices_1.push_back(ns);
            }
            std::vector<Simplex> new_simplices_2;
            for(Simplex s : new_simplices_1) {
                if (s.getVolume() < vol_eps)
                    continue;
                for(Simplex ns : s.intersectWithHyperplane(kv.first, kv.second[1])[0])
                    new_simplices_2.push_back(ns);
            }
            simplices = new_simplices_2;
        }

        double vol_prop = 0.0;
        for(Simplex s : simplices)
            vol_prop += s.getVolume() / orig_vol;
        return vol_prop;
    }
};