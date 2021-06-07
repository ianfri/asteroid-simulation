#ifndef AACELLFACE_H
#define AACELLFACE_H

#include "aaparticle.h"

class CellFace {
public:
    std::vector<std::vector<double>> vertices; // rectangle
    std::vector<double> center;
    std::vector<double> normal;
    int id;
    int nb = -1; // id of the neighboring face which has the same vertices

    bool checkIfFaceCrossPossible(Particle p, double dt);
    void calculateCenter();
};

void CellFace::calculateCenter() {
    using namespace aa;
    center = cmul(0.25, add(vertices[3], add(vertices[2], add(vertices[0], vertices[1]))));
}

bool CellFace::checkIfFaceCrossPossible(Particle p, double dt) {
    /* Uses basic ray-tracing to calculate whether or not a particle's trajectory 
    could have it cross over a face within a single timestep */

    using namespace aa;
    typedef std::vector<double> vd;

    double denominator = dot(normal, p.v);
    // std::cout << "denominator " << denominator << std::endl;
    if (std::abs(denominator) > 0.0001) { // some epsilon value to avoid working with almost orthogonal intersections
        vd difference = sub(center, p.r);
        // std::cout << "center " << center[0] << " " << center[1] << " " << center[2] << std::endl;
        double t = dot(difference, normal) / denominator;
        // std::cout << "t " << t << std::endl;
        if (t < dt) {
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

#endif