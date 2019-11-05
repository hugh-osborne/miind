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

    // std::vector<Cell> generate_cells(std::vector<unsigned int> cell_coord, std::vector<unsigned int> res) {
    //     unsigned int res_head = res[0];
        
    //     std::vector<Cell> cells;
    //     if (res.size() == 2) {
    //         std::vector<unsigned int> res_tail(res.size()-1);

    //         for (unsigned int i=0; i<res_tail.size(); i++)
    //             res_tail[i] = res[i+1];

    //         for (unsigned int d=0; d<res_head; d++) {
    //             std::vector<unsigned int> full_coord(cell_coord.size()+1);
    //             for (unsigned int c=0; c<cell_coord.size(); c++)
    //                 full_coord[c] = cell_coord[c];
    //             full_coord[cell_coord.size()] = d;

    //             std::vector<Point> points;
    //             std::vector<double> base_point_coords(num_dimensions);
    //             for (unsigned int j=0; j<num_dimensions; j++) {
    //                 base_point_coords[j] = base[j] + (full_coord[j]*(dimensions[j]/resolution[j]));
    //             }

    //             for (unsigned int i=0; i<pow(2,num_dimensions); i++) {
    //                 std::vector<double> np(num_dimensions);
    //                 for (unsigned int j=0; j<num_dimensions; j++){
    //                     if ((i >> j) & 1 == 1)
    //                         np[j] += 
    //                 }
    //                 Point p = Point()
    //                 points.push_back();
    //             }
                
    //         }
    //     }
    // }
};
