// nvcc -ccbin g++-7 -Xcompiler -O2 -Xcompiler -fopenmp -Xcompiler -std=c++14 -std=c++14 -o main

// nvcc -ccbin g++-7 -Xcompiler -O2 -Xcompiler -Wall -Xcompiler -Werror -Xcompiler -fopenmp -Xcompiler -std=c++14 -std=c++14 -o main main.cpp

// cl /EHsc /W3 /std:c++latest /openmp /O2 main.cpp

// g++-9 -O3 -march=native -fopenmp -pthread main.cpp -std=c++17 -o mainLinux -ltbb

// g++-9 -Ofast -march=native -fopenmp -pthread main.cpp -std=c++17 -o mainLinuxPar -ltbb

/* This is where the high-level code lives */
#include "aavector.h"
#include "aapopulation.h"
#include "aadomain.h"

// constexpr double INTEGRATOR_BOUND = 1.0;
constexpr double TIMELIMIT = 50;
constexpr double TIMESTEP = 0.5;
constexpr int NUM_OUTPUT_FILES = 1000; // change this when loading from backup!!
// constexpr int NUM_CELL_EDGES = 2; // actual # of cells is NUM_CELL_EDGES^3

#define PRINT false
#define LOAD_BACKUP false
const std::string LOAD_TIME = "60.500000";

// make timesteps bigger when there are no collisions, and get rid of printing when there's no collisions

int main() {

    // for now, assume we aren't reading from a backup file
        // or asking for any parameters - this is all hard-coded for now
    //double timelimit = (365.0/4.0) * 24.0 * 60.0 * 60.0; // in seconds
    // double timelimit = 15.0;
    std::cout << "[main] timelimit: " << TIMELIMIT << std::endl;

    Population population;
    population.DT = TIMESTEP;
    OutData bareData; // minimal output data
    
    double totalT;
    TimeAndParticles ntp;
    double nextTime;
    std::vector<double> timeVals;

    #if LOAD_BACKUP 
    bareData.loadBackup(population, LOAD_TIME);
    totalT = std::stod(LOAD_TIME);
    std::cout << "[main] Number of particles: " << population.pop.size() << std::endl;
    // set the index property for each particle
    for (int u = 0; u < population.pop.size(); u++) {
        population.pop[u].IDX = u;
    }
    
    ntp = population.getNextCollision_map();
    nextTime = ntp.time;
    
    double increment = (TIMELIMIT - std::stod(LOAD_TIME)) / NUM_OUTPUT_FILES;
    for (int i=1; i <= NUM_OUTPUT_FILES; i++) {
        timeVals.push_back(increment * (double)i);
    }
    #else 
    // initialize the population

    // two asteroids __ meters apart in the x-direction
    population.addAsteroid(2e11, 0.5e3, 1.03e2, {0,0,0}, {1,0,0}, 0);
    population.addAsteroid(2e11, 0.5e3, 1.03e2, {1.01e3,0,0}, {-1,0,0}, 1);

    // satellite orbit
    // population.pop.push_back(
    //     Particle({0,0,0}, {0,0,0}, 5.98e24, 6.37e6, 0)
    // );
    // population.pop.push_back(
    //     Particle({6.37e6 + 1.0e5,0,0}, {0,7.85e3,0}, 1.0, 1.0, 1)
    // );

    // create a line of particles
    // double pos = 0;
    // double vel = 0.1;
    // int color = 0;
    // for(int i=0; i < 300; i++) {
    //     population.pop.push_back(
    //         Particle({pos,0,0}, {0,vel,0}, 0.1, 0.1, color)
    //     );
    //     pos += 1;
    //     color += 1;
    // }

    // single moving asteroid
    // population.addAsteroid(2e11, 0.5e3, 1.03e2, {1.01e3,0,0}, {-1,0,0}, 0);

    // single sitting asteroid
    // population.addAsteroid(2e11, 0.5e3, 1.03e2, {0,0,0}, {0,0,0}, 0);

    // two sitting particles
    // population.pop.push_back(
    //     Particle({0,0,0}, {0,0,0}, 2e11, 0.01, 0)
    // );
    // population.pop.push_back(
    //     Particle({10,0,0}, {0,0,0}, 2e11, 0.01, 1)
    // );

    totalT = 0.0;
    bareData.addFileStreamlined(population, totalT);

    std::cout << "[main] Number of particles: " << population.pop.size() << std::endl;
    // set the index property for each particle
    for (int u = 0; u < population.pop.size(); u++) {
        population.pop[u].IDX = u;
    }
    
    ntp = population.getNextCollision_map();
    nextTime = ntp.time;
    
    double increment = TIMELIMIT / NUM_OUTPUT_FILES;
    for (int i=1; i <= NUM_OUTPUT_FILES; i++) {
        timeVals.push_back(increment * (double)i);
    }
    #endif

    int timeCounter = 1;

    std::cout << "[main] starting time: " << totalT << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(3)); 
    while (totalT < TIMELIMIT) {
        if (ntp.noColl || (!ntp.noColl && nextTime > population.DT)) { // if no collision, run one timestep
            #if PRINT
            std::cout << "[main] no collision" << std::endl;
            #endif
            population.updateAllWith4thOrderYoshida_foreach();
            
            totalT += population.DT;
            if (totalT > timeVals[timeCounter]) {
                bareData.addFileStreamlined(population, totalT);
                bareData.writeBackup(population, totalT);
                timeCounter += 1;
            }
           
            ntp = population.getNextCollision_map();
            nextTime = ntp.time;
        } else { // single collision
            #if PRINT
            std::cout << "[main] collision" << std::endl;
            #endif
            population.updateAllPositions(ntp); // handles collision as well
            
            totalT += ntp.time;
            if (totalT > timeVals[timeCounter]) {
                bareData.addFileStreamlined(population, totalT);
                bareData.writeBackup(population, totalT);
                timeCounter += 1;
            }

            ntp = population.getNextCollision_map();
            nextTime = ntp.time;
        }

        
        #if PRINT
        std::cout << "time elapsed: " << totalT << std::endl;
        std::cout << "time till next collision: " << nextTime << std::endl;
        #endif
    }
    
    bareData.addFileStreamlined(population, totalT);
    bareData.writeBackup(population, totalT);
    bareData.writeAndFlushBufferStreamlined();
}