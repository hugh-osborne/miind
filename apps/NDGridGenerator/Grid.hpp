#include "Point.hpp"
#include "Simplex.hpp"
#include "Cell.hpp"
#include "Triangulator.hpp"

class Grid {
public:
    double timestep;
    unsigned int num_dimensions;
    double threshold_v;
    Triangulator triangulator;
    std::vector<double> dimensions;
    std::vector<unsigned int> resolution;
    std::vector<double> base;
    std::vector<Point> points;
    std::vector<Point> points_trans;
    std::vector<Cell> cells;
    std::vector<Cell> cells_trans;

    Grid(std::vector<double> _base, std::vector<double> _dims, std::vector<unsigned int> _res, double _threshold_v, double _timestep):
    base(_base),
    dimensions(_dims),
    resolution(_res),
    threshold_v(_threshold_v),
    timestep(_timestep) {
        num_dimensions = _dims.size();

        // points = generate_points();
    }

    // std::vector<Point> generate_points() {

    // }
};
