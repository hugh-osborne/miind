#include "Point.hpp"
#include "Simplex.hpp"
#include "Cell.hpp"
#include "Triangulator.hpp"

#include <iostream>
#include <fstream>
#include <omp.h>
#include <iomanip>

class Grid {
public:
    double timestep;
    unsigned int num_dimensions;
    double threshold_v;
    double reset_v;
    Triangulator triangulator;
    std::vector<double> dimensions;
    std::vector<unsigned int> resolution;
    std::vector<double> base;
    std::vector<std::vector<unsigned int>> coord_list;

    Grid(std::vector<double> _base, std::vector<double> _dims, std::vector<unsigned int> _res, double _threshold_v, double _reset_v, double _timestep):
    base(_base),
    dimensions(_dims),
    resolution(_res),
    threshold_v(_threshold_v),
    reset_v(_reset_v),
    timestep(_timestep) {
        num_dimensions = _dims.size();

        generate_cell_coords(std::vector<unsigned int>(), resolution);
    }

    // obviously this constrains the number of dimensions to a hard coded 3.
    // need to pass the function as a parameter to the class. Can't be bothered right now.
    void applyRybakInterneuronEuler(Point& p) {
        double g_nap = 0.25;
        double g_na = 30.0;
        double g_k = 6.0;
        double E_na = 55.0;
        double E_k = -80.0;
        double C = 1.0;
        double g_l = 0.1;
        double E_l = -64.0;
        double I = 3.4;
        double I_h = 0.0;

        double v = p.coords[2];
        double h_na = p.coords[1];
        double m_k = p.coords[0];

        for(unsigned int i=0; i<11; i++) {

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

            v = v + (timestep/11.0)*v_prime;
            h_na = h_na + (timestep/11.0)*h_na_prime;
            m_k = m_k + (timestep/11.0)*m_k_prime;
        }

        p.coords[2] = v;
        p.coords[1] = h_na;
        p.coords[0] = m_k;

    }

    void applyHindmarshRoseEuler(Point& p) {
        double a = 1.0;
        double b = 3.0;
        double c = 1.0;
        double d = 5.0;
        double r = 0.002;
        double s = 4.0;
        double x_R = -1.6;
        double I = 3.14;
        double I_h = 0.0;

        double x = p.coords[2];
        double y = p.coords[1];
        double z = p.coords[0];

        for(unsigned int i=0; i<11; i++) {

            double x_prime = y + (-a*pow(x,3) + b*pow(x,2)) - z + I;
            double y_prime = c - (d*pow(x,2)) - y;
            double z_prime = r*(s*(x - x_R) - z);
            
            x = x + (timestep/11.0)*x_prime;
            y = y + (timestep/11.0)*y_prime;
            z = z + (timestep/11.0)*z_prime;
        }

        p.coords[2] = x;
        p.coords[1] = y;
        p.coords[0] = z;

    }

    Cell generate_cell_with_coords(std::vector<unsigned int> cell_coord, bool btranslated) {
        std::vector<double> base_point_coords(num_dimensions);
        for (unsigned int j=0; j<num_dimensions; j++) {
            base_point_coords[j] = base[j] + (cell_coord[j]*(dimensions[j]/resolution[j]));
        }

        std::vector<Point> ps = triangulator.generateUnitCubePoints(num_dimensions);
        for(unsigned int i=0; i<ps.size(); i++){
            for(unsigned int d=0; d<num_dimensions; d++){
                ps[i].coords[d] *= (dimensions[d]/resolution[d]);
                ps[i].coords[d] += base_point_coords[d];
            }
            if(btranslated)
                applyHindmarshRoseEuler(ps[i]);
        }

        return Cell(cell_coord, num_dimensions, ps, triangulator);
    }

    void generate_cell_coords(std::vector<unsigned int> cell_coord, std::vector<unsigned int> res) {
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
                full_coord[num_dimensions-1] = d;

                coord_list.push_back(full_coord);
            }

            return;
        }

        for (unsigned int d=0; d<res_head; d++){
            std::vector<unsigned int> full_coord(cell_coord.size()+1);
            for (unsigned int c=0; c<cell_coord.size(); c++)
                full_coord[c] = cell_coord[c];
            full_coord[full_coord.size()-1] = d;
            generate_cell_coords(full_coord, res_tail);
        }
    }

    std::vector<unsigned int> coords_to_strip_and_cell(std::vector<unsigned int> coords) {
        unsigned int index = 0;
        unsigned int multiplicand = 1;

        // MIIND expects membrane potential or the variable which receives instantaneous synaptic potentials
        // to be the last coordinate. So we're working back to front here.
        std::vector<unsigned int> coords_rev(num_dimensions);
        std::vector<unsigned int> res_rev(num_dimensions);
        for(unsigned int i=0; i < num_dimensions; i++) {
            // coords_rev[num_dimensions-1-i] = coords[i];
            // res_rev[num_dimensions-1-i] = resolution[i];
            coords_rev[i] = coords[i];
            res_rev[i] = resolution[i];
        }

        for (unsigned int res : res_rev) 
            multiplicand *= res;
        for(int d=0; d<num_dimensions-1; d++) {
            multiplicand /= res_rev[d];
            index += int(coords_rev[d] * multiplicand);
        }
        index /= res_rev[num_dimensions-1];
        std::vector<unsigned int> pair(2);
        pair[0] = index;
        pair[1] = coords_rev[num_dimensions-1];
        return pair;
    }

    unsigned int coords_to_index(std::vector<unsigned int> coords) {
        unsigned int index = 0;
        unsigned int operand = 1;
        for (int c=num_dimensions-1; c>=0; c--) {
            index += coords[c] * operand;
            operand *= resolution[c];
        }
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

    void buildCellRange(std::vector<Cell>& cell_ptrs, std::vector<unsigned int> base_min, std::vector<unsigned int> max_coords, std::vector<unsigned int> min_coords) {

        if (max_coords.size() == 1) {
            for(unsigned int i=0; i<(max_coords[0] - min_coords[0])+1; i++) {
                std::vector<unsigned int> nb = base_min;
                nb.push_back(min_coords[0] + i);
                cell_ptrs.push_back(generate_cell_with_coords(nb, false));
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

    std::vector<Cell> getCellRange(Cell& tcell) {

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

        std::vector<Cell> cells;
        buildCellRange(cells, std::vector<unsigned int>(), max_coords, min_coords);
        return cells;
    }
    
    std::map<std::vector<unsigned int>,double>
    calculateTransitionForCell(Cell& tcell, std::vector<Cell>& cell_range) {
        std::map<std::vector<unsigned int>,double> t;
        for(Cell check_cell : cell_range) {
            double prop = tcell.intersectsWith(check_cell);
            if (prop == 0)
                continue;
            t[check_cell.grid_coords] = prop;
        } 
        return t;
    }

    std::map<std::vector<unsigned int> ,std::map<std::vector<unsigned int>, double>> calculateTransitionMatrix() {  
        std::map<std::vector<unsigned int> ,std::map<std::vector<unsigned int>, double>> transitions;    
#pragma omp parallel for
        for (unsigned int c=0; c < coord_list.size(); c++) {
            Cell cell = generate_cell_with_coords(coord_list[c],true);
            std::vector<Cell> check_cells = getCellRange(cell);
            std::map<std::vector<unsigned int>, double> ts = calculateTransitionForCell(cell, check_cells);

            if (ts.size() == 0) { // cell was completely outside the grid, so don't move it.
                ts[cell.grid_coords] = 1.0;
            }

            double total_prop = 0.0;
            for(auto const& kv : ts){
                total_prop += kv.second;
            }
            double missed_prop = 1.0/total_prop;

            for(auto const& kv : ts){
                double d = ts[kv.first];
                ts[kv.first] *= missed_prop;
            }
            transitions[cell.grid_coords] = ts;
            if(transitions.size() % 100 == 0)
                std::cout << transitions.size() << " complete(ish).\n";
        }
        return transitions;
    }

    void generateTMatFileBatched(std::string basename) { 
        unsigned int batch_size = 1000;
        std::ofstream file;
        file.open(basename + ".tmat");

        file << "0\t0\n";

        std::map<std::vector<unsigned int> ,std::map<std::vector<unsigned int>, double>> transitions;  
        for (unsigned int batch=0; batch < coord_list.size() / batch_size; batch++) {
#pragma omp parallel for
            for (unsigned int c=(batch*batch_size); c < (batch*batch_size)+batch_size; c++) {
                Cell cell = generate_cell_with_coords(coord_list[c],true);
                std::vector<Cell> check_cells = getCellRange(cell);
                std::map<std::vector<unsigned int>, double> ts = calculateTransitionForCell(cell, check_cells);

                if (ts.size() == 0) { // cell was completely outside the grid, so don't move it.
                    ts[cell.grid_coords] = 1.0;
                }

                double total_prop = 0.0;
                for(auto const& kv : ts){
                    total_prop += kv.second;
                }
                double missed_prop = 1.0/total_prop;

                for(auto const& kv : ts){
                    double d = ts[kv.first];
                    ts[kv.first] *= missed_prop;
                }
#pragma omp critical
                transitions[cell.grid_coords] = ts;
            }

            for(auto const& kv : transitions) {
                std::vector<unsigned int> pair = coords_to_strip_and_cell(kv.first);
                file << "1000000000;" << pair[0] << "," << pair[1] << ";";
                for(auto const& tv : kv.second) {
                    std::vector<unsigned int> tpair = coords_to_strip_and_cell(tv.first);
                    file << tpair[0] << "," << tpair[1] << ":" << tv.second << ";";
                }
                file << "\n";
            }

            std::cout << '\r' << std::setw(5) << 100.0 * ((float)(batch*batch_size)/(float)coord_list.size()) << "% complete." << std::setfill(' ') << std::flush;
            transitions.clear();
        }
        std::cout << "\n";
    }

    void generateTMatFileLowMemory(std::string basename) {
        std::ofstream file;
        file.open(basename + ".tmat");

        file << "0\t0\n";
        unsigned int num_lines = 0;
        for(std::vector<unsigned int> c : coord_list) {
            Cell cell = generate_cell_with_coords(c,true);
            std::vector<unsigned int> pair = coords_to_strip_and_cell(cell.grid_coords);
            file << "1000000000;" << pair[0] << "," << pair[1] << ";";
            std::vector<Cell> check_cells = getCellRange(cell);
            std::map<std::vector<unsigned int>, double> ts = calculateTransitionForCell(cell, check_cells);

            if (ts.size() == 0) { // cell was completely outside the grid, so don't move it.
                ts[cell.grid_coords] = 1.0;
            }

            double total_prop = 0.0;
            for(auto const& kv : ts){
                total_prop += kv.second;
            }
            double missed_prop = 1.0/total_prop;

            for(auto const& kv : ts){
                double d = ts[kv.first];
                ts[kv.first] *= missed_prop;
                std::vector<unsigned int> tpair = coords_to_strip_and_cell(kv.first);
                file << tpair[0] << "," << tpair[1] << ":" << kv.second << ";";
            }
            
            file << "\n";
            num_lines++;
            if(num_lines % 100 == 0)
                std::cout << num_lines << " complete.\n";
        }

        file.close();
    }

    void generateTMatFile(std::string basename) {

        std::map<std::vector<unsigned int> ,std::map<std::vector<unsigned int>, double>> trs = calculateTransitionMatrix();

        std::ofstream file;
        file.open(basename + ".tmat");

        file << "0\t0\n";
        for(auto const& kv : trs) {
            std::vector<unsigned int> pair = coords_to_strip_and_cell(kv.first);
            file << "1000000000;" << pair[0] << "," << pair[1] << ";";
            for(auto const& tv : kv.second) {
                std::vector<unsigned int> tpair = coords_to_strip_and_cell(tv.first);
                file << tpair[0] << "," << tpair[1] << ":" << tv.second << ";";
            }
            file << "\n";
        }

        file.close();
    }

    void generateModelFile(std::string basename) {
        std::ofstream file;
        file.open(basename + ".model");

        file << "<Model>\n";
        file << "<Mesh>\n";
        file << "<TimeStep>" << timestep*0.001 << "</TimeStep>\n";
        file << "</Mesh>\n";
        file << "<Stationary>\n";
        file << "</Stationary>\n";
        file << "<Mapping type = \"Reversal\">\n";
        file << "</Mapping>\n";
        file << "<threshold>" << threshold_v << "</threshold>\n";
        file << "<V_reset>" << reset_v << "</V_reset>\n";

        file << "<Mapping type=\"Reset\">\n";
        file << "</Mapping>\n";

        file << "</Model>\n";

        file.close();
    }

};
