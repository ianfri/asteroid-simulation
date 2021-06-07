#ifndef AATIMEPART_H
#define AATIMEPART_H

#include "aavector.h"
#include "aaparticle.h" // for general includes

struct TimeAndParticles{
        double time = 0.0; // time of collision between two particles
        std::vector<int> particleIndices = {0,0};

        bool noColl = true;

        bool operator() (TimeAndParticles a, TimeAndParticles b) { // for std::sort
            return a.time < b.time;
        }
};

#endif