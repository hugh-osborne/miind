#ifndef APP_ND_GRID_POINT
#define APP_ND_GRID_POINT

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

    Point() :
    coords(0),
    lives(0),
    connected(0),
    dead(false),
    hyper(false) {}

    Point& operator=(const Point &other) {
        coords = std::vector<double>(other.coords.size());
        for(unsigned int i=0; i<other.coords.size(); i++)
            coords[i] = other.coords[i];

        connected = std::vector<Point*>(other.connected.size());
        for(unsigned int i=0; i<other.connected.size(); i++)
            connected[i] = other.connected[i];

        lives = other.lives;
        dead = other.dead;
        hyper = other.hyper;

        return *this;
    }

    bool operator==(const Point &other) const {
        if(other.coords.size() != coords.size())
            return false;
        
        bool equal = true;
        for(unsigned int i=0; i<coords.size(); i++){
            equal &= (coords[i] == other.coords[i]);
        }
        return equal;
    }

};

#endif