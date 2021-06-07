#ifndef AACELL_H
#define AACELL_H

#include "aacellface.h"
#include "aaparticle.h"

class Cell {
public:
    std::vector<int> indices; // indices of the particles presently in the cell
    std::vector<CellFace> faces; // rectangular prism

    // EventQueue eQ;

    int checkIfCellCrossPossible(Particle p, double dt); // -1 means no cross, otherwise face id
    void solve();
    // void solve(std::vector<Particle> pList);
};

int Cell::checkIfCellCrossPossible(Particle p, double dt) {
    for (auto face : faces) {
        if (face.checkIfFaceCrossPossible(p, dt)) {
            return face.id;
        }
    }
    return -1;
}

#endif