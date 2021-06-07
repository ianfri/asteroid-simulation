#ifndef AAPARTICLE_H
#define AAPARTICLE_H

// general includes that propagate upwards
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip> // setprecision
#include <cmath>
#include <random>
#include <limits>
#include <cfenv>
#include <execution>
#include <iostream>
#include <set>
#include <thread>

#include "aavector.h"

class Particle {
// Data
public:
    std::vector<double> r; // position components
    std::vector<double> v; // velocity components
    double m;
    double radius;
    int color;
    int cellId;
    bool crossOver = false;
    int IDX; // index; doesn't change

public:
    // Constructor
    Particle(std::vector<double> r_, std::vector<double> v_, double mass_, double radius_, int color_) {
        r = r_;
        v = v_;
        m = mass_;
        radius = radius_;
        color = color_;
    }

    // Methods
    double calculateCollisionTime(Particle p); // return time of collision with Particle p
    void collide(Particle &p); // calculate new velocities; MODIFIES INPUT PARAM

    // Particle operator= (Particle other) {
    //     for (int i=0; i < this->r.size(); i++) {
    //         if 
    //     }
    // }
};

double Particle::calculateCollisionTime(Particle p) {
    double b = aa::dot( aa::sub(p.r, r), aa::sub(p.v, v) );
    double u = aa::mag( aa::sub(p.v, v) );
    double dist = aa::mag( aa::sub(r, p.r) );
    double d = 2.0*radius;
    double discr = b*b - u*u*(dist*dist - d*d);
    return (-1.0*b - sqrt(discr)) / (u*u);

    // test for doubleing point error
    if ( std::fetestexcept(FE_ALL_EXCEPT) ) {
        std::cout << "<<WARNING>> [calculateCollisionTime] doubleing-point error" << std::endl;
        std::feclearexcept(FE_ALL_EXCEPT);
    }
}

void Particle::collide(Particle &p) {
    // note: assumes same radius and mass for both
    double b = aa::dot( aa::sub(p.r, r), aa::sub(p.v, v) ); // know is good
    double d = 2.0*radius;
    std::vector<double> delV = aa::cmul(1.0 / (d*d),  aa::cmul(b, aa::sub(p.r, r)));
    v = aa::add(v, delV);
    p.v = aa::sub(p.v, delV);
}

#endif
