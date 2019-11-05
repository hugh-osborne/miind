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
    std::vector<Cell> cells;
    std::vector<Cell> cells_trans;

    Grid(std::vector<double> _base, std::vector<double> _dims, std::vector<unsigned int> _res, double _threshold_v, double _timestep):
    base(_base),
    dimensions(_dims),
    resolution(_res),
    threshold_v(_threshold_v),
    timestep(_timestep) {
        num_dimensions = _dims.size();

        generate_cells(std::vector<unsigned int>(), resolution);
    }

    // obviously this constrains the number of dimensions to a hard coded 3.
    // need to pass the function as a parameter to the class. Can't be bothered right now.
    void applyFunctionEuler(Point& p) {
        double g_nap = 0.25;
        double g_na = 30.0;
        double g_k = 6.0;
        double E_na = 55.0;
        double E_k = -80.0;
        double C = 1.0;
        double g_l = 0.1;
        double E_l = -64.0;
        double I = 0.15;
        double I_h = 0.0;

        double v = p.coords[2];
        double h_na = p.coords[1];
        double m_k = p.coords[0];

        double I_l = -g_l*(v - E_l);
        double I_na = -g_na * h_na * (v - E_na) * (pow((pow((1 + exp(-(v+35)/7.8)),-1)),3));
        double I_k = -g_k * (pow(m_k,4)) * (v - E_k);

        double v_prime = I_na + I_k + I_l + I;

        double part_1 = (pow((1 + (exp((v + 55)/7))),-1)) - h_na;
        double part_2 = 30 / (exp((v+50)/15) + exp(-(v+50)/16));
        double h_na_prime = part_1 / part_2;

        double part_1 = (pow((1 + (exp(-(v + 28)/15))),-1)) - m_k;
        double part_2 = 7 / (exp((v+40)/40) + exp((-v + 40)/50));
        double m_k_prime = part_1 / part_2;

        p.coords[0] = m_k + (timestep)*m_k_prime;
        p.coords[1] = h_na + (timestep)*h_na_prime;
        p.coords[2] = v + (timestep)*v_prime;

    }

    void generate_cells(std::vector<unsigned int> cell_coord, std::vector<unsigned int> res) {
        unsigned int res_head = res[0];

        if (res.size() <= 1) // We don't work with 1D
            return;

        std::vector<unsigned int> res_tail(res.size()-1);

        for (unsigned int i=0; i<res_tail.size(); i++)
            res_tail[i] = res[i+1];
        
        if (res.size() == 2) {
            for (unsigned int d=0; d<res_head; d++) {
                std::vector<unsigned int> full_coord(cell_coord.size()+1);
                for (unsigned int c=0; c<cell_coord.size(); c++)
                    full_coord[c] = cell_coord[c];
                full_coord[cell_coord.size()] = d;

                std::vector<double> base_point_coords(num_dimensions);
                for (unsigned int j=0; j<num_dimensions; j++) {
                    base_point_coords[j] = base[j] + (full_coord[j]*(dimensions[j]/resolution[j]));
                }

                std::vector<Point> points(pow(2,num_dimensions));
                for (unsigned int i=0; i<pow(2,num_dimensions); i++) {
                    std::vector<double> np(num_dimensions);
                    for (unsigned int j=0; j<num_dimensions; j++){
                        if ((i >> j) & 1 == 1)
                            np[j] = base_point_coords[j] + (dimensions[j]/resolution[j]);
                        else
                            np[j] = base_point_coords[j];
                    }
                    Point p = Point(np);
                    points.push_back(p);
                }
                cells.push_back(Cell(full_coord, num_dimensions, points, triangulator));
            }

            return;
        }

        for (unsigned int d=0; d<res_head; d++){
            std::vector<unsigned int> full_coord(cell_coord.size()+1);
            for (unsigned int c=0; c<cell_coord.size(); c++)
                full_coord[c] = cell_coord[c];
            full_coord[cell_coord.size()] = d;
            generate_cells(full_coord, res_tail);
        }
    }
    
};
