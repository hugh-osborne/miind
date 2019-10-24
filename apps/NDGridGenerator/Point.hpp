#include <vector>

class Point {
public:
    std::vector<double> coords;
    std::vector<Point*> connected;
    int lives;
    bool dead;
    bool hyper;

    Point(std::vector<double> _coords) :
    coords(_coords),
    lives(0),
    connected(0),
    dead(false),
    hyper(false) {}

    bool Point::operator==(const Point &other) const {
        if(other.coords.size() != coords.size())
            return false;
        
        bool equal = true;
        for(unsigned int i=0; i<coords.size(); i++){
            equal &= (coords[i] == other.coords[i]);
        }
        return equal;
    }

};