#ifndef AADOMAIN_H
#define AADOMAIN_H

// testing out the function for creating 6 different planes

// g++-9 -I SUPPORT/boost_UNIX -O3 -march=native -fopenmp -pthread testing_boxes.cpp -std=c++17 -o testing_boxes -ltbb

#include <vector>
#include <iostream>
#include <set>
#include <thread>

#include "aacell.h"
#include "aapopulation.h"
#include "aaoutdata.h"

// #include <boost/fusion/algorithm/transformation/flatten.hpp>
// #include <boost/fusion/include/flatten.hpp>
// #include <boost/fusion/container/vector/convert.hpp>
// #include <boost/fusion/include/as_vector.hpp>

// for use in for_each
void Cell_solve(Cell &cell) {
    std::cout << "hello" << std::endl;
    for (auto i : cell.indices) std::cout << i << " ";
    std::cout << std::endl;
}

class Domain {
public:
    std::vector<Cell> cells;
    std::vector< std::vector<std::pair<int, int>> > crossingIds; // {particle id, id of face it could cross}
    // EventQueue eQCross; // event queue for the crossovers

    void create(Population pop, int numCells); // for now, numCells isn't actually the number of cells
    void assignNeighboringFaces();
    void checkCrossovers(Population pop, double dt); // modify crossingIds
    void solve(Population pop, OutData &dataFile); // uses crossingIds
};

// void Domain::checkCrossovers(Population pop, double dt) {
//     for (auto cell : cells) {
//         int pCounter = 0;
//         for (auto p : pop.pop) {
//             int tmp = cell.checkIfCellCrossPossible(p, dt);
//             if (tmp != -1) {
//                 crossingIds.push_back({pCounter, tmp});
//             }
//             pCounter += 1;
//         }
//     }
// }

// change this to aggreate the results for each cell sepeartely instead of using
    // the global vector crossingIds
void Domain::checkCrossovers(Population pop, double dt) {
    // TODO: crossingIds = (something below)
    // std::for_each (cells.begin(), cells.end(),  [&](Cell cell) {
    // std::vector<std::pair<int, int>> a(cells.size());
    std::vector< std::vector<std::pair<int, int>> > a (cells.size());
    std::transform (cells.begin(), cells.end(), a.begin(),  [&](Cell cell) {
        int pCounter = 0;
        std::vector<std::pair<int, int>> tmpCrossingIds;
        for (auto p : pop.pop) {
            int tmp = cell.checkIfCellCrossPossible(p, dt);
            if (tmp != -1) {
                // CHANGE
                tmpCrossingIds.push_back({pCounter, tmp});
            }
            pCounter += 1;
        }
        return tmpCrossingIds;
    });

    // should be small, so copying is okay?
    crossingIds = a;

    // crossingIds = boost::fusion::as_vector(boost::fusion::flatten(a));
}

void Domain::create(Population pop, int numCells) {
    std::vector<double> mins = {0, 0, 0}; // x,y,z
    std::vector<double> maxs = {5, 5, 5};
    for (int i=0; i<3; i++) { mins[i] -= 1; maxs[i] += 1; }

    // calculate the amount that will be added to each coordinate each time a new point is added to the domain
    int numDivisions = numCells; // TODO: fix in the future, derive numDivisions from numCells
    std::vector<double> offsets;
    for (int i=0; i<3; i++) { 
        double tmp = (std::abs(mins[i]) + std::abs(maxs[i])) / (double)numDivisions;
        std::cout << tmp << std::endl;
        offsets.push_back(tmp);
    }

    // REMEMBER THAT THESE MULTIPLIERS ARE FOR THE BOTTOM SMALLEST CORNER OF EACH BOX

    std::vector<std::vector<int>> multipliers;
    for (int xm=0; xm < numDivisions; xm++) {
        for (int ym=0; ym < numDivisions; ym++) {
            for (int zm=0; zm < numDivisions; zm++) {
                multipliers.push_back({xm, ym, zm});
            }
        }
    }
    
    std::vector<double> origMins = {mins[0], mins[1], mins[2]};
    for (auto mult : multipliers) {

        /* Create a temporary cell and then add it to cells[] */
        std::vector<CellFace> tmpFaces(6);
        // xy
        tmpFaces[0].vertices.push_back(
            {mins[0], mins[1], mins[2]}
        );
        tmpFaces[0].vertices.push_back(
            {mins[0] + offsets[0], mins[1], mins[2]}
        );
        tmpFaces[0].vertices.push_back(
            {mins[0], mins[1] + offsets[1], mins[2]}
        );
        tmpFaces[0].vertices.push_back(
            {mins[0] + offsets[0], mins[1] + offsets[1], mins[2]}
        );
        tmpFaces[0].normal = {0,0,1};

        // xy + z
        tmpFaces[1].vertices.push_back(
            {mins[0], mins[1], mins[2] + offsets[2]}
        );
        tmpFaces[1].vertices.push_back(
            {mins[0] + offsets[0], mins[1], mins[2] + offsets[2]}
        );
        tmpFaces[1].vertices.push_back(
            {mins[0], mins[1] + offsets[1], mins[2] + offsets[2]}
        );
        tmpFaces[1].vertices.push_back(
            {mins[0] + offsets[0], mins[1] + offsets[1], mins[2] + offsets[2]}
        );
        tmpFaces[1].normal = {0,0,-1};

        // yz
        tmpFaces[2].vertices.push_back(
            {mins[0], mins[1], mins[2]}
        );
        tmpFaces[2].vertices.push_back(
            {mins[0], mins[1], mins[2] + offsets[2]}
        );
        tmpFaces[2].vertices.push_back(
            {mins[0], mins[1] + offsets[1], mins[2]}
        );
        tmpFaces[2].vertices.push_back(
            {mins[0], mins[1] + offsets[1], mins[2] + offsets[2]}
        );
        tmpFaces[2].normal = {1,0,0};

        // yz + x
        tmpFaces[3].vertices.push_back(
            {mins[0] + offsets[0], mins[1], mins[2]}
        );
        tmpFaces[3].vertices.push_back(
            {mins[0] + offsets[0], mins[1], mins[2] + offsets[2]}
        );
        tmpFaces[3].vertices.push_back(
            {mins[0] + offsets[0], mins[1] + offsets[1], mins[2]}
        );
        tmpFaces[3].vertices.push_back(
            {mins[0] + offsets[0], mins[1] + offsets[1], mins[2] + offsets[2]}
        );
        tmpFaces[3].normal = {-1,0,0};

        // xz
        tmpFaces[4].vertices.push_back(
            {mins[0], mins[1], mins[2]}
        );
        tmpFaces[4].vertices.push_back(
            {mins[0], mins[1], mins[2] + offsets[2]}
        );
        tmpFaces[4].vertices.push_back(
            {mins[0] + offsets[0], mins[1], mins[2]}
        );
        tmpFaces[4].vertices.push_back(
            {mins[0] + offsets[0], mins[1], mins[2] + offsets[2]}
        );
        tmpFaces[4].normal = {0,1,0};

        // xz + y
        tmpFaces[5].vertices.push_back(
            {mins[0], mins[1] + offsets[1], mins[2]}
        );
        tmpFaces[5].vertices.push_back(
            {mins[0], mins[1] + offsets[1], mins[2] + offsets[2]}
        );
        tmpFaces[5].vertices.push_back(
            {mins[0] + offsets[0], mins[1] + offsets[1], mins[2]}
        );
        tmpFaces[5].vertices.push_back(
            {mins[0] + offsets[0], mins[1] + offsets[1], mins[2] + offsets[2]}
        );
        tmpFaces[5].normal = {0,-1,0};

        int idx = 0;
        for (auto i : tmpFaces) {
            std::cout << "FACE NO. " << idx << std::endl;
            for (auto j : i.vertices) {
                for (auto k : j) {
                    std::cout << k << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
            idx += 1;
        }

        Cell tmpCell;
        for (auto face : tmpFaces) {
            tmpCell.faces.push_back(face);
        }
        for (auto &face: tmpCell.faces) {
            face.calculateCenter();
        }
        cells.push_back(tmpCell);

        mins[0] = origMins[0] + mult[0]*offsets[0];
        mins[1] = origMins[1] + mult[1]*offsets[1];
        mins[2] = origMins[2] + mult[2]*offsets[2];

    }
    if (multipliers.size() != cells.size()) {
        std::cout << "<<WARNING>>[Domain::create] multipliers.size() != cells.size()" << std::endl;
    }
    std::cout << "<<INFO>>[Domain::create] number of cells: " << cells.size() << std::endl;

    /* Now, load the population into each of the cells */

    for (auto &cell : cells) {
        int particleCounter = 0;
        for (auto p : pop.pop) {
            // "if the x component of p is greater than the x component of the lower left vertice of the bottom face" etc...
            if (p.r[0] > cell.faces[0].vertices[0][0] && p.r[0] < cell.faces[0].vertices[1][0]) {
                if (p.r[1] > cell.faces[0].vertices[0][1] && p.r[1] < cell.faces[0].vertices[2][1]) {
                    if (p.r[2] > cell.faces[4].vertices[0][2] && p.r[2] < cell.faces[4].vertices[1][2]) {
                        cell.indices.push_back(particleCounter);
                    }
                }
            }
        }
        particleCounter += 1;
    }

    int sanity = 0;
    for (auto cell : cells) {
        sanity += (int)cell.indices.size();
    }
    std::cout << "<<INFO>>[Domain::create] number of particles in cells: " << sanity << std::endl;
           
}

void Domain::assignNeighboringFaces() {

    // assign id's to each face
    int idx = 0;
    for (auto &cell : cells) {
        for (auto &face : cell.faces) { // do normal indexing INSTEAD
            face.id = idx;
            idx += 1;
        }
    }

    auto pairIfNeighbors = [&](CellFace a, CellFace b) {
        std::set aPoints = {a.vertices[0], a.vertices[1], a.vertices[2], a.vertices[3]};
        std::set bPoints = {b.vertices[0], b.vertices[1], b.vertices[2], b.vertices[3]};
        if (aPoints == bPoints && a.id != b.id) {
            std::cout << "neighbor " << a.id << " " << b.id << std::endl;
            a.nb = b.id;
            b.nb = a.id;
        }
    };

    std::vector<CellFace> allFacesTmp;
    for (auto cell : cells) {
        for (auto face : cell.faces) {
            allFacesTmp.push_back(face);
        }
    }

    for (int i=0; i < allFacesTmp.size(); i++) {
        for (int j=0; j < allFacesTmp.size(); j++) {
            pairIfNeighbors(allFacesTmp[i], allFacesTmp[j]); 
        }
    }

    idx = 0;
    for (auto &cell : cells) {
        for (auto &face : cell.faces) {
            face.nb = allFacesTmp[idx].nb;
            idx += 1;
        }
    }
}

void Domain::solve(Population pop, OutData &dataFile) {

    // check for crossovers

    // PASS IN THE INFO ABOUT CROSSOVERS

    checkCrossovers(pop, pop.DT);

    std::for_each(cells.begin(), cells.end(), Cell_solve);

    // std::vector<std::thread> threadPool;
    // if (cells.size() < numThreads) {
    //     // execute on less than numThreads
    // } else if (cells.size() % numThreads == 0) {
    //     // split up the work equally
    //     int cellsPerThread = cells.size() % numThreads;
        

    // } else {
    //     int leftover = cells.size() % numThreads;

    // }

}

void mainDomain() {

    Population popu;

    popu.pop.push_back(
        Particle({1,1,1}, {1.0, 0.1, 0.1}, 5, 1, 0)
    );
    popu.pop.push_back(
        Particle({5,5,5}, {0,0,0}, 5, 1, 1)
    );

    Domain domain;
    domain.create(popu, 2);

    std::cout << "Face cross possible?: " <<
    domain.cells[0].faces[3].checkIfFaceCrossPossible(popu.pop[0], 1.4) << std::endl;

    domain.assignNeighboringFaces();

    OutData theOut;
    domain.solve(popu, theOut);

}

#endif