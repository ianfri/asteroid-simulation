#ifndef AAPOPULATION_H
#define AAPOPULATION_H

#include "aacell.h"
#include "aaparticle.h"
#include "aatimepart.h"

#include <execution>

#include <mutex>
std::mutex protect0;
std::mutex protect1;

struct AsteroidCoordinate {
    std::vector<double> r;
    bool valid; // flag for if valid coordinate
};

class Population {
public:
// Data
    std::vector<Particle> pop;
    double DT; // standard timestep
    double G = 	6.67430e-11;//6.673e-11;//6.67430e-11; //std::pow(6.674, -11.0);
    int numCells;
    // std::vector<Cell> cells;

// Methods
    // add a homogenous random cube of particles to the population for benchmarking
    void addRandPop(int popSize, int bound, double maxVel);
    // add an asteroid to the population
        // prad is particle radius
    void addAsteroid(double massAsteroid, double radius, 
                     double prad, std::vector<double> r_, std::vector<double> v_, int color);
    TimeAndParticles getNextCollision();
    TimeAndParticles getNextCollision_map();
    TimeAndParticles getBestCollisionForParticle(Particle theParticle);
    // void calculateAllCollisionTimes_new(EventQueue &eventQ); // use cells and smarter event queue
    // void calculateAllCollisionTimesCell (Cell &cell);
    // void calculateAllCollisionTimesCellIncoming (Cell &cell, int incomingIdx);
     // move each particle on a straight path until the next collision in the queue at time time
    void updateAllPositions(TimeAndParticles nextTimeAndParticles);
    // integrators; choose one
    void updateAllWithRK4(); // NOT working; particles barely move at all
    void updateAllWithVVerlet(); // 2nd-order symplectic
    void updateAllWith4thOrderYoshida(); // 4th-order symplectic
    void updateAllWith4thOrderYoshida_foreach();
    std::vector<double> calcAcceleration(Particle &p, std::vector<double> currentR, double totalM);

// Grid methods
    // void partitionDomain();
    // void putParticlesInCells();
};

void Population::addRandPop(int popSize, int bound, double maxVel) {
    if (popSize > (bound*bound*bound)) {
        std::cout << "<<ERROR>>[Population::addRandPop] popsize > volume" << std::endl;
    }
    
    std::random_device rd;     // only used once to initialize (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(0, bound);

    std::random_device rdf;     // only used once to initialize (seed) engine
    std::mt19937 rngf(rdf());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_real_distribution<double> unif( -maxVel, maxVel );

    // create the population
    for (int i = 0; i < popSize; i++) {
        std::vector<double> rID = { (double)uni(rng), (double)uni(rng), (double)uni(rng) };
        std::vector<double> vID = { unif(rngf), unif(rngf), unif(rngf) };

        // mass of 1kg, radius of 0.5kg
        Particle tempP(rID, vID, 1.0, 0.5, 0);
        pop.push_back(tempP);
    }
}

void Population::addAsteroid(double massAsteroid, double radius, double prad, std::vector<double> r_, std::vector<double> v_, int color) {
    double buffer = 0.001; // extra distance between particles in the same asteroid
    double rRange = 1.5*radius; //generous
    std::vector<double> rMins = aa::csub(rRange, r_);
    std::vector<double> rMaxs = aa::cadd(rRange, r_);

    // populate the asteroid grid, start in the bottom left-hand corner
    // TODO: please make this one (or two?) loops by using the fact that rMins & rMaxs are vectors!
    std::vector<AsteroidCoordinate> asteroidGrid;
    for (double currentX = rMins[0]; currentX < rMaxs[0]; currentX = currentX + prad + buffer) {
        for (double currentY = rMins[1]; currentY < rMaxs[1]; currentY = currentY + prad + buffer) {
            for (double currentZ = rMins[2]; currentZ < rMaxs[2]; currentZ = currentZ + prad + buffer) {
                AsteroidCoordinate temp;
                temp.r = {currentX, currentY, currentZ};
                temp.valid = true;
                asteroidGrid.push_back(temp);
            }
        }
    }

    // find the valid coordinates
    int numInvalidPoints = 0;
    for (unsigned int i=0; i < asteroidGrid.size(); i++) {
        std::vector<double> diff = aa::sub(asteroidGrid.at(i).r, r_);
        double coordinateDistance = aa::mag(diff);
        if ( coordinateDistance > radius ) {
            asteroidGrid.at(i).valid = false;
            numInvalidPoints++;
        }
    }
    int numValidPoints = (int)asteroidGrid.size() - numInvalidPoints;

    // make the particles in the asteroid
    double pMass = massAsteroid / numValidPoints;
    for (unsigned int i=0; i < asteroidGrid.size(); i++) {
        if ( asteroidGrid.at(i).valid ) {
            Particle tempP(asteroidGrid.at(i).r, v_, pMass, prad, color);
            pop.push_back(tempP);
        }
    }
}

// TO PARALLELIZE
    // probably won't be able to use Thrust easily
    // but could use CPU parallel
// TODO: don't calculate duplicate collisions!! ~2x speedup
// TODO: this is first priority for parallelization.

/* Think about using std::transform

There could be a function, let's call it getBestCollisionForParticle(Particle theParticle)

std::transform(pop.begin(), pop.end(), outputList.begin(), getBestCollisionForParticle());

then just loop through outputList to find the best one and return it
*/

// TimeAndParticles Population::getBestCollisionForParticle(Particle theParticle) {
//     TimeAndParticles best;
//     best.time = std::numeric_limits<double>::max();
//     int i = theParticle.IDX;
//     for (unsigned int j=0; j < pop.size(); j++) {
//         double dst = aa::mag( aa::sub( pop[i].r, pop[j].r ) );  
//         double filter_0 = (pop[i].radius + pop[j].radius); 
//         double filter_1 = 1.1 * DT * (aa::mag(pop[i].v) + aa::mag(pop[j].v));
//         double filter = filter_0 + filter_1;

//         if ( (i != j) && (dst < filter) ) {
//             double tempTime = pop[i].calculateCollisionTime( pop[j] );
//             if ((tempTime > 0) && (tempTime < best.time)) {
//                 best.particleIndices = {i,j};
//                 best.time = tempTime;
//                 best.noColl = false;
//             }
//         } 
//     }
//     return best;
// }

TimeAndParticles Population::getNextCollision_map() {
    TimeAndParticles best;
    best.time = std::numeric_limits<double>::max();
    std::vector<TimeAndParticles> outputList(pop.size());
    std::transform(std::execution::par, pop.begin(), pop.end(), outputList.begin(), [this] (Particle theParticle) {
        TimeAndParticles bestInner;
        bestInner.time = std::numeric_limits<double>::max();
        int i = theParticle.IDX;
        for (int j=0; j < pop.size(); j++) {
            double dst = aa::mag( aa::sub( pop[i].r, pop[j].r ) );  
            double filter_0 = (pop[i].radius + pop[j].radius); 
            double filter_1 = 1.1 * DT * (aa::mag(pop[i].v) + aa::mag(pop[j].v));
            double filter = filter_0 + filter_1;

            if ( (i != j) && (dst < filter) ) {
                double tempTime = pop[i].calculateCollisionTime( pop[j] );
                if ((tempTime > 0) && (tempTime < bestInner.time)) {
                    bestInner.particleIndices = {i,j};
                    bestInner.time = tempTime;
                    bestInner.noColl = false;
                }
            } 
        }
        return bestInner;
    });
    for (auto tp : outputList) {
        if (tp.time < best.time && !tp.noColl) {
            best.time = tp.time;
            best.particleIndices = tp.particleIndices;
            best.noColl = false;
        }
    }

    return best;
}

TimeAndParticles Population::getNextCollision() {
    TimeAndParticles best;
    best.time = std::numeric_limits<double>::max();
    for (int i=0; i < pop.size(); i++) {
        for (int j=0; j < pop.size(); j++) {
            // CHANGE TO COMPARE IF i!=j FIRST

            // ALSO, try using the distance squared to compare for the filter instead
                // to eliminate the square root!
                // consider replacing aa::mag with aa:mag2 (for magnitude squared)
            double dst = aa::mag( aa::sub( pop[i].r, pop[j].r ) );  
            double filter_0 = (pop[i].radius + pop[j].radius); 
            double filter_1 = 1.1 * DT * (aa::mag(pop[i].v) + aa::mag(pop[j].v));
            double filter = filter_0 + filter_1;

            if ( (i != j) && (dst < filter) ) {
                double tempTime = pop[i].calculateCollisionTime( pop[j] );
                if ((tempTime > 0) && (tempTime < best.time)) {
                    best.particleIndices = {i,j};
                    best.time = tempTime;
                    best.noColl = false;
                }
            } 
           
        }
    }
    return best;
}

/*
// just within the cell
void Population::calculateAllCollisionTimesCell (Cell &cell) {
    int eventCounter = 0;

    // only make the OUTER for_each parallel
    std::for_each(cell.indices.begin(), cell.indices.end(), [this, &cell, &eventCounter](int i) {
        std::for_each(cell.indices.begin(), cell.indices.end(), [this, &cell, i, &eventCounter](int j) {
            double tempTime;
            if (i == j) {
                tempTime = std::numeric_limits<double>::max() / 2.0;
            } else {  
                double dst = aa::mag( aa::sub( pop[i].r, pop[j].r ) );  
                // double filter = (pop[i].radius + pop[j].radius) + 1.1 * DT * (aa::mag(pop[i].v) + aa::mag(pop[j].v));
                double filter_0 = pop[i].radius + pop[j].radius; 
                double filter_1 = 1.1 * DT * (aa::mag(pop[i].v) + aa::mag(pop[j].v)); // 1.1 for extra safety
                double filter = filter_0 + filter_1;

                if ( dst < filter ) {
                    tempTime = pop[i].calculateCollisionTime( pop[j] );
                } else {
                    tempTime = std::numeric_limits<double>::max() / 2.0;
                }
                if (eventCounter < cell.eQ.queue.size()) {
                    cell.eQ.queue[eventCounter].time = tempTime;
                    cell.eQ.queue[eventCounter].particleIndices[0] = i;
                    cell.eQ.queue[eventCounter].particleIndices[1] = j;
                } else {
                    TimeAndParticles tmpTaP;
                    tmpTaP.time = tempTime;
                    tmpTaP.particleIndices[0] = i;
                    tmpTap.particleIndices[1] = j;
                    cell.eQ.queue.push_back(tmpTaP);
                }
                eventCounter += 1;
            }
        });

        // if there are leftover queue elements from the previous timestep cancel them; shortening vectors is expensive
        if (cell.eQ.queue.size() / eventCounter != 1) {
            std::for_each(cell.eQ.queue.end() - (cell.eQ.queue.size() % eventCounter), cell.eQ.queue.end(), [&](TimeAndParticles tap) {
                tap.time = std::numeric_limits<double>::max() / 2.0;
            });
        }
    });
}

// a singular outside particle with all the particles within the cell
// CALL AFTER!!! calculateAllCollisionTimesCell
void Population::calculateAllCollisionTimesCellIncoming (Cell &cell, int incomingIdx) {
    int i = incomingIdx;
    std::for_each(cell.indices.begin(), cell.indices.end(), [this, i](int j) {
        if (i == j) {
            tempTime = std::numeric_limits<double>::max() / 2.0;
        } else {
            double tempTime;
            double dst = aa::mag( aa::sub( pop[i].r, pop[j].r ) );  
            // double filter = (pop[i].radius + pop[j].radius) + 1.1 * DT * (aa::mag(pop[i].v) + aa::mag(pop[j].v));
            double filter_0 = pop[i].radius + pop[j].radius; 
            double filter_1 = 1.1 * DT * (aa::mag(pop[i].v) + aa::mag(pop[j].v)); // 1.1 for extra safety
            double filter = filter_0 + filter_1;

            if ( dst < filter ) {
                tempTime = pop[i].calculateCollisionTime( pop[j] );
            } else {
                tempTime = std::numeric_limits<double>::max() / 2.0;
            }
            
            // avoid allocating new memory if possible; look for a dud TimeAndParticles to overwrite
            bool overwrote = false;
            for (auto &tap : cell.eQ.queue) {
                if (!overwrote && std::abs((std::numeric_limits<double>::max() / 2.0) - tap.time) < 0.01) {
                    tap.time = tmpTime;
                    tap.particleIndices[0] = i;
                    tap.particleIndices[1] = j;
                    overwrote = true;
                }
            }
            if (!overwrote) {
                TimeAndParticles tmpTaP;
                tmpTaP.time = tempTime;
                tmpTaP.particleIndices[0] = i;
                tmpTap.particleIndices[1] = j;
                cell.eQ.queue.push_back(tmpTaP);
            }
        }
    });
}
*/

/*
void Population::calculateAllCollisionTimes_new(EventQueue &eventQ) {
    int eventCounter = 0;

    auto collsInCell = [this](Particle p, Cell c, EventQueue &eQ, int &eventCounter) {
        for (auto p_ : c.indices) {
            double tempTime;
            double dst = aa::mag( aa::sub( p.r, pop[p_].r ) );  
            
            double filter_0 = (p.radius + pop[p_].radius); 
            double filter_1 = 1.1 * DT * (aa::mag(p.v) + aa::mag(pop[p_].v));

            double filter = filter_0 + filter_1;

            if ( (p.IDX != p_) && (dst < filter) ) {
                tempTime = p.calculateCollisionTime( pop[p_] );
            } else {
                tempTime = std::numeric_limits<double>::max() / 2.0;
            }

            if (eventCounter < eQ.queue.size()) { // queue allocated enough
                eQ.queue[eventCounter].time = tempTime;
                eQ.queue[eventCounter].particleIndices[0] = p.IDX;
                eQ.queue[eventCounter].particleIndices[1] = p_; 
                eventCounter += 1;
            } else { // queue not allocated enough
                TimeAndParticles tmp;
                eQ.queue.push_back(tmp);
                eQ.queue[eventCounter].time = tempTime;
                eQ.queue[eventCounter].particleIndices[0] = p.IDX;
                eQ.queue[eventCounter].particleIndices[1] = p_; 
                eventCounter += 1;
            }
        }
    };

    for (auto &part : pop) {
        int crossCellId = cells[part.cellId].checkIfCellCrossPossible(part);
        if (crossCellId == -1) {
            // just check within cell
            collsInCell(part, cells[part.cellId], eventQ, eventCounter);
        } else if (crossCellId >= 0 && crossCellId < cells.size()) {
            // check within cell and the neighboring crossover cell
            collsInCell(part, cells[part.cellId], eventQ, eventCounter);
            collsInCell(part, cells[crossCellId], eventQ, eventCounter);
        } else {
            // the particle must have moved OUTSIDE the grid
                // no collisions, but still do gravity
        }
    }

    // clean up any extras in the event queue
    if (eventCounter < eventQ.queue.size()) {
        for (int i = eventCounter+1; i < eventQ.queue.size(); i++) {
            eventQ.queue[i].time = std::numeric_limits<double>::max() / 2.0;
        }
    }

}
*/

void Population::updateAllPositions(TimeAndParticles nextTimeAndParticles) {
    double nextTime = nextTimeAndParticles.time;
    int pID_0 = nextTimeAndParticles.particleIndices.at(0);
    int pID_1 = nextTimeAndParticles.particleIndices.at(1);

    // First, move the particles up to the next collision time
    // TO PARALLELIZE
        // could easily use transform() or foreach()
    for (unsigned int i=0; i < pop.size(); i++) {
        pop[i].r = aa::add(pop[i].r, aa::cmul(nextTime, pop[i].v));
    }

    // Next, deal with the collision
    pop[pID_0].collide( pop[pID_1] );
}


// ADD MUTEXES HERE?????
std::vector<double> Population::calcAcceleration(Particle &p, std::vector<double> currentR, double totalM) {
    typedef std::vector<double> vd;
    using namespace aa;
    vd centerR = {0,0,0};
    for (Particle i : pop) {
        if (p.IDX != i.IDX) {
            vd a = cmul(i.m / totalM, i.r);
            centerR = add(centerR, a);
        }
    }
    std::vector<double> r_cm_p = aa::sub(centerR, currentR); // vector p -> center of mass
    return aa::cmul(
            (G * totalM) / mag2(r_cm_p), 
            aa::cmul(1.0 / mag(r_cm_p), r_cm_p)
            );
}

void Population::updateAllWith4thOrderYoshida_foreach() {
    // for convenience
    typedef std::vector<double> vd;
    using namespace aa;

    // make sure that the average velocity is 0
    vd avgVel = {0, 0, 0};
    for (auto p : pop) {
        avgVel = add(avgVel, cmul(1.0/pop.size(), p.v));
    }
    std::cout << "AVG SPEED: " << mag(avgVel) << std::endl;

    // for (Particle &p : pop) {
    std::for_each(std::execution::par, pop.begin(), pop.end(), [this] (Particle &p) {
        // calculate the center of mass of all of the particles except the current one
        double totalM = 0.0;
        for (Particle i : pop) {
            if (p.IDX != i.IDX) {
                if (totalM + i.m < totalM) {
                    std::cout << "<<WARNING>> [updateAllWith4thOrderYoshida] totalM may have tried to overflow" << std::endl;
                } else if (std::isnan(totalM + i.m)) {
                    std::cout << "<<WARNING>> [updateAllWith4thOrderYoshida] totalM tried to be NaN" << std::endl;
                } else {
                    totalM += i.m;
                }
            }
        }

        vd x1, x2, x3, v1, v2, v3;
        const double w_0 = -1.0 * pow(2.0, 1.0 / 3.0) / (2.0 - pow(2.0, 1.0 / 3.0));
        const double w_1 = 1.0 / (2.0 - pow(2.0, 1.0 / 3.0));
        const double c_1 = w_1 / 2.0;
        const double c_4 = c_1;
        const double c_2 = (w_0 + w_1) / 2.0;
        const double c_3 = c_2;
        const double d_1 = w_1;
        const double d_3 = d_1;
        const double d_2 = w_0;

        std::vector<double> centerR = {0,0,0};
        // protect0.lock();
        for (Particle i : pop) {
            if (p.IDX != i.IDX) {
                vd tmp = cmul(i.m / totalM, i.r);
                centerR = add(centerR, tmp);
            }
        }
        // protect0.unlock();
        
        // calculate gravitational acceleration a
        std::vector<double> r_cm_p = aa::sub(centerR, p.r); // vector p -> center of mass
        
        x1 = add(p.r, cmul(c_1 * DT, p.v));
        v1 = add(p.v, cmul(d_1 * DT, this->calcAcceleration(p, x1, totalM)));
        // std::cout << "[yosh] firstAcc: " << mag(cmul(d_1 * DT, this->calcAcceleration(p, x1, totalM))) << std::endl;
        x2 = add(x1, cmul(c_2 * DT, v1));
        v2 = add(v1, cmul(d_2 * DT, this->calcAcceleration(p, x2, totalM)));
        x3 = add(x2, cmul(c_3 * DT, v2));
        v3 = add(v2, cmul(d_3 * DT, this->calcAcceleration(p, x3, totalM)));

        vd tmpR = add(x3, cmul(c_4 * DT, v3));
        vd tmpV = v3;

        // protect1.lock();
        // p.r = add(x3, cmul(c_4 * DT, v3));
        // p.v = v3; 
        p.r = tmpR;
        p.v = tmpV;
        // protect1.unlock();
    });
    // }
}

void Population::updateAllWith4thOrderYoshida() {
    // for convenience
    typedef std::vector<double> vd;
    using namespace aa;

    for (Particle &p : pop) {
        // calculate the center of mass of all of the particles except the current one
        double totalM = 0.0;
        for (Particle i : pop) {
            if (p.IDX != i.IDX) {
                if (totalM + i.m < totalM) {
                    std::cout << "<<WARNING>> [updateAllWith4thOrderYoshida] totalM may have tried to overflow" << std::endl;
                } else if (std::isnan(totalM + i.m)) {
                    std::cout << "<<WARNING>> [updateAllWith4thOrderYoshida] totalM tried to be NaN" << std::endl;
                } else {
                    totalM += i.m;
                }
            }
        }
        std::vector<double> centerR = {0,0,0};
        for (Particle i : pop) {
            if (p.IDX != i.IDX) {
                vd tmp = cmul(i.m / totalM, i.r);
                centerR = add(centerR, tmp);
            }
        }


        // calculate gravitational acceleration a
        std::vector<double> r_cm_p = aa::sub(centerR, p.r); // vector p -> center of mass
        const double w_0 = -1.0 * pow(2.0, 1.0 / 3.0) / (2.0 - pow(2.0, 1.0 / 3.0));
        const double w_1 = 1.0 / (2.0 - pow(2.0, 1.0 / 3.0));
        const double c_1 = w_1 / 2.0;
        const double c_4 = c_1;
        const double c_2 = (w_0 + w_1) / 2.0;
        const double c_3 = c_2;
        const double d_1 = w_1;
        const double d_3 = d_1;
        const double d_2 = w_0;

        vd x1, x2, x3, v1, v2, v3;
        
        x1 = add(p.r, cmul(c_1 * DT, p.v));
        v1 = add(p.v, cmul(d_1 * DT, this->calcAcceleration(p, x1, totalM)));
        x2 = add(x1, cmul(c_2 * DT, v1));
        v2 = add(v1, cmul(d_2 * DT, this->calcAcceleration(p, x2, totalM)));
        x3 = add(x2, cmul(c_3 * DT, v2));
        v3 = add(v2, cmul(d_3 * DT, this->calcAcceleration(p, x3, totalM)));
        p.r = add(x3, cmul(c_4 * DT, v3));
        p.v = v3; 

    }
}

void Population::updateAllWithVVerlet() {
    // for convenience
    typedef std::vector<double> vd;
    using namespace aa;

    // calculate the center of mass of all of the particles
    // only do this once for each integrator call
        // okay assumption for now; the center of mass will barely change with each integrator call
    double totalM = 0.0;
    for (Particle i : pop) {
        if (totalM + i.m < totalM) {
            std::cout << "<<WARNING>> [updateAllWithVVerlet] totalM may have tried to overflow" << std::endl;
        } else if (std::isnan(totalM + i.m)) {
            std::cout << "<<WARNING>> [updateAllWithVVerlet] totalM tried to be NaN" << std::endl;
        } else {
            totalM += i.m;
        }
    }

    std::vector<double> centerR = {0,0,0};
    for (Particle i : pop) {
        volatile double proportion = i.m / totalM;
        centerR = aa::cmul(proportion, aa::add(i.r, centerR));
        if (std::isnan(centerR[0]) || std::isnan(centerR[1]) || std::isnan(centerR[2])) {
            std::cout << "<<WARNING>> [updateAllWithVVerlet] a component of centerR is NaN" << std::endl;
        } 
    }

    // Vectorized RK4
    for (Particle &p : pop) {
    //concurrency::parallel_for_each( pop.begin(), pop.end(), [&](Particle p) {
        // calculate gravitational acceleration a
        std::vector<double> r_cm_p = aa::sub(centerR, p.r); // vector p -> center of mass
        double mag_r_cmp_p = aa::mag(r_cm_p);
        std::vector<double> a = aa::cmul(
            (G * totalM) / (mag_r_cmp_p * mag_r_cmp_p), 
            aa::cmul(1.0 / mag_r_cmp_p, r_cm_p)
            );
        
       // Velocity Verlet
       // https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet

        vd nPos = add(p.r, add(cmul(DT, p.v), cmul(DT*DT*0.5, a)));
        vd nVel = add(p.v, cmul(2.0*DT*0.5, a)); // a bit of a cheat, not using new_acc b/c barely changes

        p.r = nPos;
        p.v = nVel;
    }
}

// TO PARALLELIZE:
    // turn into std::transform or foreach after
void Population::updateAllWithRK4() {
    // for convenience
    typedef std::vector<double> vd;
    using namespace aa;

    // calculate the center of mass of all of the particles
    // only do this once for each RK4 call
        // okay assumption for now that the center of mass will barely change with each RK4 call
    double totalM = 0.0;
    for (Particle i : pop) {
        if (totalM + i.m < totalM) {
            std::cout << "<<WARNING>> [updateAllWithRK4] totalM may have tried to overflow" << std::endl;
        } else if (std::isnan(totalM + i.m)) {
            std::cout << "<<WARNING>> [updateAllWithRK4] totalM tried to be NaN" << std::endl;
        } else {
            totalM += i.m;
        }
    }

    std::vector<double> centerR = {0,0,0};
    // for (Particle i : pop) {
    //     volatile double proportion = i.m / totalM;
    //     centerR = aa::cmul(proportion, aa::add(i.r, centerR));
    //     // centerR = aa::add(centerR, aa::cmul(proportion, i.r));
    //     if (std::isnan(centerR[0]) || std::isnan(centerR[1]) || std::isnan(centerR[2])) {
    //         std::cout << "<<WARNING>> [updateAllWithRK4] a component of centerR is NaN" << std::endl;
    //     } 
    //     //totalM += i.m;
    // }
    for (Particle i : pop) {
        centerR = add(centerR, cmul(i.m, i.r));
    }
    centerR = cmul(1.0 / totalM, centerR);
    //centerR = aa::cmul( 1.0 / totalM, centerR );

    // std::cout << "CENTER OF MASS\n" << centerR[0] << " " << centerR[1] << " " << centerR[2] << std::endl;
    // std::cout << "TOTAL MASS\n" << totalM << std::endl;

    // Vectorized RK4
    for (Particle &p : pop) {
    //concurrency::parallel_for_each( pop.begin(), pop.end(), [&](Particle p) {
        // calculate gravitational acceleration a
        std::vector<double> r_cm_p = aa::sub(centerR, p.r); // vector p -> center of mass
        // std::vector<double> r_cm_p = aa::sub(p.r, centerR);
        //std::cout << "checkpoint 1" << std::endl;
        double mag_r_cmp_p = aa::mag(r_cm_p);
        std::vector<double> a = aa::cmul(
            (G * totalM) / (mag_r_cmp_p * mag_r_cmp_p), 
            aa::cmul(1.0 / mag_r_cmp_p, r_cm_p)
            );

        // std::cout << "ACCELERATION\n" << a[0] << " " << a[1] << " " << a[2] << std::endl;
        
        // calculate k-coefficients
        // NOTE: WE ARE LETTING H = 1.0
        vd k1_r = p.v;
        vd k1_v = compprod(a, p.r);

        vd k2_r = cmul(DT*0.5, compprod(p.v, k1_v));
        vd k2_v = compprod(a, add(p.r, cmul(DT*0.5, k1_r)));

        vd k3_r = cmul(DT*0.5, compprod(p.v, k2_v));
        vd k3_v = compprod(a, add(p.r, cmul(DT*0.5, k2_r)));

        vd k4_r = compprod(p.v, cmul(DT, k3_v));
        vd k4_v = compprod(a, add(p.r, cmul(DT, k3_r)));

        // update velocity and position
        // CHECK THIS: adding in delta t right at the end... is this okay??
        vd kv_sum = add(k1_v, add(k4_v, add(cmul(2.0, k2_v), cmul(2.0, k3_v))));
        // std::cout << "VELOCITY\n"<< "p.v: " << p.v[0] << " " << p.v[1] << " " << p.v[2] << std::endl;
        p.v = add(p.v, cmul((1.0/6.0) * DT, kv_sum));
        // std::cout << "p.v: " << p.v[0] << " " << p.v[1] << " " << p.v[2] << std::endl;
        // std::cout << "kv_sum: " << kv_sum[0] << " " << kv_sum[1] << " " << kv_sum[2] << std::endl;

        vd kr_sum = add(k1_r, add(k4_r, add(cmul(2.0, k2_r), cmul(2.0, k3_r))));
        // std::cout << "POSITION\n" << "p.r: " << p.r[0] << " " << p.r[1] << " " << p.r[2] << std::endl;
        p.r = add(p.r, cmul((1.0/6.0) * DT, kr_sum));
        // std::cout << "p.r: " << p.r[0] << " " << p.r[1] << " " << p.r[2] << std::endl;
        // std::cout << "kr_sum: " << kr_sum[0] << " " << kr_sum[1] << " " << kr_sum[2] << std::endl;
    //});
    }
}

// void Population::partitionDomain() {
//     // get the maximum maginitude x, y, and z coordinates of the pop
//     std::vector<double> mags = {0, 0, 0};
//     for (auto p : pop) {
//         for (int i = 0; i < p.r.size(); i++) {
//             if (abs(p.r[i]) > mags[i]) {
//                 mags[i] = abs(p.r[i]);
//             }
//         }
//     }
//     for (auto i : mags) { std::cout << "<<DEBUG>> mag " << i << std::endl; }

//     std::vector<double> edgeLengths = aa::cmul(2.0 / numCells, mags); // 2 b/c -'s; giant cube
//     for (auto i : edgeLengths) { std::cout << "<<DEBUG>> edge " << i << std::endl; }
    
//     double xTally = -mags[0];
//     double yTally = -mags[1];
//     double zTally = -mags[2];

//     std::vector<std::vector<double>> tmpPoints;

//     std::vector<double> xVals;
//     std::vector<double> yVals;
//     std::vector<double> zVals;


//     // FIND THE X-Y-Z values that specify the lower left corners of the domain
//     for (int count = 0; count < numCells+1; count++) { 
//         Cell cell;
//         CellFace face;
//     }



//     for (int count = 0; count < numCells+1; count++) {
//         xVals.push_back(xTally + (count * edgeLengths[0]));
//         yVals.push_back(xTally + (count * edgeLengths[1]));
//         zVals.push_back(xTally + (count * edgeLengths[2]));
//     }


//     for (int xCount = 0; xCount < numCells+1; xCount++) {
//         for (int yCount = 0; yCount < numCells+1; yCount++) {
//             for (int zCount = 0; zCount < numCells+1; zCount++) {
//                 tmpPoints.push_back({
//                     xTally + (xCount * edgeLengths[0]),
//                     yTally + (yCount * edgeLengths[1]),
//                     zTally + (zCount * edgeLengths[2])
//                 });
//             }
//         }
//     }

//     for (auto i : tmpPoints) {
//         std::cout << "<<DEBUG>> point ";
//         for (auto j : i) {
//             std::cout << j << " ";
//         }
//         std::cout << std::endl;
//     }

// }

#endif