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

        generate_cells(std::vector<unsigned int>(), resolution, false);
        generate_cells(std::vector<unsigned int>(), resolution, true);
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

        double v = p.coords[0];
        double h_na = p.coords[1];
        double m_k = p.coords[2];

        double I_l = -g_l*(v - E_l);
        double I_na = -g_na * h_na * (v - E_na) * (pow((pow((1 + exp(-(v+35)/7.8)),-1)),3));
        double I_k = -g_k * (pow(m_k,4)) * (v - E_k);

        double v_prime = I_na + I_k + I_l + I;

        double part_1 = (pow((1 + (exp((v + 55)/7))),-1)) - h_na;
        double part_2 = 30 / (exp((v+50)/15) + exp(-(v+50)/16));
        double h_na_prime = part_1 / part_2;

        part_1 = (pow((1 + (exp(-(v + 28)/15))),-1)) - m_k;
        part_2 = 7 / (exp((v+40)/40) + exp((-v + 40)/50));
        double m_k_prime = part_1 / part_2;

        p.coords[2] = m_k + (timestep)*m_k_prime;
        p.coords[1] = h_na + (timestep)*h_na_prime;
        p.coords[0] = v + (timestep)*v_prime;

    }

    void generate_cells(std::vector<unsigned int> cell_coord, std::vector<unsigned int> res, bool btranslated) {
        unsigned int res_head = res[0];

        if (res.size() < 1) // We don't work with 1D
            return;

        std::vector<unsigned int> res_tail(res.size()-1);

        for (unsigned int i=0; i<res_tail.size(); i++)
            res_tail[i] = res[i+1];

        if (res.size() == 1) {
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
                    if(btranslated)
                        applyFunctionEuler(p);
                    points[i] = p;
                }
                if(btranslated)
                    cells_trans.push_back(Cell(full_coord, num_dimensions, points, triangulator));
                else
                    cells.push_back(Cell(full_coord, num_dimensions, points, triangulator));
            }

            return;
        }

        for (unsigned int d=0; d<res_head; d++){
            std::vector<unsigned int> full_coord(cell_coord.size()+1);
            for (unsigned int c=0; c<cell_coord.size(); c++)
                full_coord[c] = cell_coord[c];
            full_coord[cell_coord.size()] = d;
            generate_cells(full_coord, res_tail, btranslated);
        }
    }

    unsigned int coords_to_index(std::vector<unsigned int> coords) {
        unsigned int index = 0;
        unsigned int operand = 1;
        for (int c=num_dimensions-1; c>=0; c--) {
            index += coords[c] * operand;
            operand *= resolution[c];
        }
        std::vector<unsigned int> co = index_to_coords(index);
        return index;
    }

    std::vector<unsigned int> index_to_coords(unsigned int index) {
        std::vector<unsigned int> coords(num_dimensions);
        unsigned int operand = 1;

        std::vector<unsigned int> operands(num_dimensions);

        for (int c=num_dimensions-1; c>=0; c--) {
            operands[c] = operand;
            operand *= resolution[c];
        }

        for (int c=0; c<num_dimensions-1; c++) {
            unsigned int m = int(index / operands[c]);
            coords[c] = m;
            index = index % operands[c];
        }
        coords[num_dimensions-1] = index;
        return coords;
    }

    std::vector<unsigned int> getCellCoordsForPoint(Point& p){
        std::vector<unsigned int> coords(num_dimensions);
        for(unsigned int c=0; c<num_dimensions; c++) {
            int co = int((p.coords[c] - base[c]) / (dimensions[c]/resolution[c]));
            if (co >= int(resolution[c]))
                co = resolution[c]-1;
            if (co < 0)
                co = 0;
            coords[c] = co;
        }
        return coords;
    }

    void buildCellRange(std::vector<Cell*>& cell_ptrs, std::vector<unsigned int> base_min, std::vector<unsigned int> max_coords, std::vector<unsigned int> min_coords) {

        if (max_coords.size() == 1) {
            for(unsigned int i=0; i<(max_coords[0] - min_coords[0])+1; i++) {
                std::vector<unsigned int> nb = base_min;
                nb.push_back(min_coords[0] + i);
                cell_ptrs.push_back(&cells[coords_to_index(nb)]);
            }
        } else  if (max_coords.size() > 1) {
            for(unsigned int i=0; i<(max_coords[0] - min_coords[0])+1; i++){
                std::vector<unsigned int> nb = base_min;
                nb.push_back(min_coords[0] + i);

                std::vector<unsigned int> max_tail(max_coords.size()-1);    
                for(unsigned int j=0; j<max_coords.size()-1; j++)
                    max_tail[j] = max_coords[j+1];

                std::vector<unsigned int> min_tail(min_coords.size()-1);    
                for(unsigned int j=0; j<min_coords.size()-1; j++)
                    min_tail[j] = min_coords[j+1];

                buildCellRange(cell_ptrs, nb, max_tail, min_tail);
            }
        }
    }

    std::vector<Cell*> getCellRange(Cell& tcell) {
        std::vector<double> max = tcell.simplices[0].points[0].coords;
        std::vector<double> min = tcell.simplices[0].points[0].coords;

        for (unsigned int c=0; c<num_dimensions; c++) {
            for (Simplex s : tcell.simplices){
                for (Point p : s.points) {
                    if (max[c] < p.coords[c])
                        max[c] = p.coords[c];
                    if (min[c] > p.coords[c])
                        min[c] = p.coords[c];
                }
            }
        }

        Point p_max = Point(max);
        Point p_min = Point(min);

        std::vector<unsigned int> max_coords = getCellCoordsForPoint(p_max);
        std::vector<unsigned int> min_coords = getCellCoordsForPoint(p_min);

        std::vector<Cell*> cells;
        buildCellRange(cells, std::vector<unsigned int>(), max_coords, min_coords);
        return cells;
    }
    
    std::map<std::vector<unsigned int>,double>
    calculateTransitionForCell(Cell& tcell, std::vector<Cell*> cell_range) {
        std::map<std::vector<unsigned int>,double> t;
        std::cout << tcell.grid_coords[0] << " " << tcell.grid_coords[1] << " " << tcell.grid_coords[2] << ": \n";
        for(Cell* check_cell : cell_range) {
            double prop = tcell.intersectsWith(*check_cell);
            if (prop == 0)
                continue;
            t[check_cell->grid_coords] = prop;
        } 
        return t;
    }

    void calculateTransitionMatrix() {      
        for (Cell cell : cells_trans) {
            std::vector<Cell*> check_cells = getCellRange(cell);
            std::map<std::vector<unsigned int>, double> ts = calculateTransitionForCell(cell, check_cells);

            if (ts.size() == 0) { // cell was completely outside the grid, so don't move it.
                ts[cell.grid_coords] = 1.0;
            }

            double total_prop = 0.0;
            for(auto const& kv : ts){
                total_prop += kv.second;
            }
            double missed_prop = 1.0 - total_prop;
            double share_prop = missed_prop / ts.size();

            for(auto const& kv : ts){
                ts[kv.first] += share_prop;
                std::cout << kv.first[0] << "," << kv.first[1] << "," << kv.first[2] << " : " << ts[kv.first] << "\n";
            }
        }
    }

};
