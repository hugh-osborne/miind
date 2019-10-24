#include "Point.hpp"
#include "Triangulator.hpp"

#include <stdlib.h>
#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/io.hpp> 
#include <boost/numeric/ublas/matrix_proxy.hpp> 
#include <boost/numeric/ublas/lu.hpp> 

class Simplex {
public:
    Triangulator& triangulator;
    unsigned int num_dimensions;
    std::vector<Point> points;
    std::vector<Point> lines;

    Simplex(unsigned int num_dims, std::vector<std::vector<double>> _points, Triangulator& _triangulator):
    num_dimensions(num_dims),
    points(_points.size()),
    lines(0),
    triangulator(_triangulator) {
        for (unsigned int i=0; i<_points.size(); i++) {
            points[i] = Point(_points[i]);
        }

        for (unsigned int i=0; i<points.size(); i++) {
            points[i].connected = std::vector<Point*>(points.size()-1);
            for (unsigned int j=0; j<points.size(); j++) {
                if (j == i) continue;
                points[i].connected[j] = &points[j];
            }
        }

        lines = generateLines();
    }

    std::vector<Point> generateLines() {
        std::vector<Point> lines(num_dimensions);
        for(unsigned int p=0; p<points.size()-1; p++) {
            std::vector<double> coords(num_dmiensions);
            for (unsigned int c=0; c<num_dimensions; c++)
                coords[c] = points[p].coords[c] - points[0].coords[0];
            lines[p] = Point(coords);
        }
        return lines;
    }

    // CalcDeterminant by Richel Bilderbeek : http://www.richelbilderbeek.nl/CppUblasMatrixExample7.htm
    double CalcDeterminant(boost::numeric::ublas::matrix<double> m) 
    { 
        assert(m.size1() == m.size2() && "Can only calculate the determinant of square matrices"); 
        boost::numeric::ublas::permutation_matrix<std::size_t> pivots(m.size1() ); 

        const int is_singular = boost::numeric::ublas::lu_factorize(m, pivots); 

        if (is_singular) return 0.0; 

        double d = 1.0; 
        const std::size_t sz = pivots.size(); 
        for (std::size_t i=0; i != sz; ++i) 
        { 
            if (pivots(i) != i) 
            { 
            d *= -1.0; 
            } 
            d *= m(i,i); 
        } 
        return d; 
    } 

    double getVolume() {
        boost::numeric::ublas::matrix<double> m(num_dimensions,num_dimensions);
        for (unsigned int l=0; l<num_dimensions; l++) {
            for(unsigned int c=0; c<num_dimensions; c++) {
                m(l,c) = lines[l].coords[c];
            }
        }

        unsigned int dim_fac = 0
        for(unsigned int n=0; n<num_dimensions; n++)
            dim_fac += n;

        return abs(CalcDeterminant(m)/dim_fac);
    }

    std::vector<std::vector<Simplex>> intersectWithHyperplane(unsigned int dim_index, double dim) {
        double eps = 0.00000000001;

        std::vector<Point*> lower;
        std::vector<Point*> upper;
        for (Point p : points) {
            if(p.coords[dim_index] < dim - eps) lower.push_back(&p);
            if(p.coords[dim_index] > dim + eps) upper.push_back(&p);
        }

        std::vector<Point> p_outs;
        for (Point* p0 : lower){
            for (Point* p1 : upper) {
                double t = (dim - p0->coords[dim_index]) / (p1->coords[dim_index] - p0->coords[dim_index]);
                std::vector<double> coords(num_dimensions);
                for (unsigned int i=0; i<num_dimensions; i++){
                    coords[i] = p0->coords[i] + ((p1->coords[i] - p0->coords[i])*t);
                }
                Point np(coords);
                std::vector<Point*> cs(2);
                cs[0] = &p0;
                cs[1] = &p1;
                np.connected = cs;
                for (Point* p : p0->connected)
                    if (*p == *p1) p = &np;
                for (Point* p : p1->connected)
                    if (*p == *p0) p = &np;
                np.hyper = true;
                p_outs.push_back(np);

            }
        }

        if (p_outs.size() == 0) {
            std::vector<std::vector<Simplex>> out(2);
            if (points[0].coords[dim_index] > dim){
                std::vector<Simplex> less(1);
                less[0] = Simplex(num_dimensions, points, triangulator);
                out[0] = less;
                out[1] = std::vector<Simplex>(0);
            } else if (points[0].coords[dim_index] < dim){
                std::vector<Simplex> greater(1);
                greater[0] = Simplex(num_dimensions, points, triangulator);
                out[0] = std::vector<Simplex>(0);
                out[1] = greater; 
            }     
            return out;
        }

        unsigned int index = 0;
        std::vector<unsigned int> i_less(lower.size());
        for (unsigned int i=0; i<lower.size(); i++) i_less[i] = i + index;
        index += lower.size();

        std::vector<unsigned int> i_greater(upper.size());
        for (unsigned int i=0; i<upper.size(); i++) i_greater[i] = i + index;
        index += upper.size();

        std::vector<unsigned int> i_hyp(p_outs.size());
        for (unsigned int i=0; i<p_outs.size(); i++) i_hyp[i] = i + index;
        index += p_outs.size();

        std::vector<Point*> p_equal;
        for (Point p : points) {
            if (p.coords[dim_index] <= dim + eps && p.coords[dim_index] >= dim - eps)
                p_equal.push_back(&p);
        }
        
        for (unsigned int i=0; i<p_equal.size(); i++) i_hyp.push_back(i + index);

        std::vector<Point*> p_total(lower.size()+upper.size()+p_outs.size()+p_equal.size());
        for(unsigned int i=0; i<lower.size(); i++){
            p_total[i] = lower[i];
        }
        for(unsigned int i=0; i<upper.size(); i++){
            p_total[lower.size()+i] = upper[i];
        }
        for(unsigned int i=0; i<p_outs.size(); i++){
            p_total[lower.size()+upper.size()+i] = &p_outs[i];
        }
        for(unsigned int i=0; i<p_equal.size(); i++){
            p_total[lower.size()+upper.size()+p_outs.size()+i] = p_equal[i];
        }

        std::vector<Simplex*> simplices = triangulator.chooseTriangulation(p_total, i_less, i_greater, i_hyp);

        std::vector<Simplex> less;
        std::vector<Simplex> greater;
        for (Simplex* s : simplices){
            bool all_above = true;
            bool all_below = true;
            for (Point p : s->points) {
                all_above &= p.coords[dim_index] >= dim;
                all_below &= p.coords[dim_index] <= dim;
            }

            if (all_above)
                greater.push_back(*s);

            if (all_below)
                less.push_back(*s);
        }

        std::vector<std::vector<Simplex>> out(2);
        out[0] = less;
        out[1] = greater;
        return out;
    }
};