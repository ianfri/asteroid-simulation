#ifndef AAOUTDATA_H
#define AAOUTDATA_H

#include "aapopulation.h"

constexpr auto OUTPUT_DIR = "OUT/TESTS/two/";

class OutData {
public:
    std::vector<std::string> buffersStreamlined;

    // all info - basically convert the Population to a string and write to disk
    void writeBackup(Population population, double timestep);

    void loadBackup(Population &population, const std::string timestep);

    // less info - basically convert the Population to a string and store in a vector
    void addFileStreamlined(Population population, double timestep);
    void writeAndFlushBufferStreamlined();

};

void OutData::writeBackup(Population population, double timestep) {
    /*
    * TIMESTEP,ID,RX,RY,RZ,VX,VY,VZ,MASS,RADIUS,COLOR
    */
   std::string tempString = "";
   std::ostringstream stringStream;
   stringStream << std::fixed;
   stringStream << std::setprecision(16);
   for (unsigned int i=0; i < population.pop.size(); i++) {
       stringStream <<  timestep;
       stringStream <<  ",";
       stringStream <<  i;
       stringStream <<  ",";
       stringStream <<  population.pop.at(i).r[0];
       stringStream <<  ",";
       stringStream <<  population.pop.at(i).r[1];
       stringStream <<  ",";
       stringStream <<  population.pop.at(i).r[2];
       stringStream <<  ",";
       stringStream <<  population.pop.at(i).v[0];
       stringStream <<  ",";
       stringStream <<  population.pop.at(i).v[1];
       stringStream <<  ",";
       stringStream <<  population.pop.at(i).v[2];
       stringStream <<  ",";
       stringStream <<  population.pop.at(i).m;
       stringStream <<  ",";
       stringStream <<  population.pop.at(i).radius;
       stringStream <<  ",";
       stringStream <<  population.pop.at(i).color;
       stringStream <<  "\n";
   }
    std::ofstream outputFile;
    outputFile.open(OUTPUT_DIR + std::to_string(1) + "backups/" + std::to_string(timestep) + "_backup.csv");
    outputFile << stringStream.str();
    outputFile.close();
}

void OutData::loadBackup(Population &population, const std::string timestep) {
    std::ifstream inputFile;
    inputFile.open(OUTPUT_DIR + std::to_string(1) + "backups/" + timestep + "_backup.csv");
    if (inputFile.is_open()) {
        std::string line;
        while (std::getline(inputFile, line)) {
            std::stringstream lineStream(line);
            std::vector<std::string> sV;
            while (lineStream.good()) {
                std::string tmp;
                std::getline(lineStream, tmp, ',');
                sV.push_back(tmp);
            }

            using namespace std;
            population.pop.push_back(
                Particle( {stod(sV[2]), stod(sV[3]), stod(sV[4])}, {stod(sV[5]), stod(sV[6]), stod(sV[7])},
                    stod(sV[8]), stod(sV[9]), stoi(sV[10])
                )
            );
        }
    } else {
        std::cout << "[OutData::loadBackup] <<ERROR>> unable to open backup" << std::endl;
    }
}

void OutData::addFileStreamlined(Population population, double timestep) {
    std::cout << "[OutData::addFileStreamlined] Added position data for time __" << timestep << "__ to buffer" << std::endl;
    /*
    * TIMESTEP,RX,RY,RZ
    */
   std::string tempString = "";
   std::ostringstream stringStream;
   stringStream << std::fixed;
   stringStream << std::setprecision(16);
   stringStream << "time," << "rx," << "ry," << "rz," << "color" << "\n";
   for (unsigned int i=0; i < population.pop.size(); i++) {
       stringStream <<  timestep;
       stringStream <<  ",";
       stringStream <<  population.pop.at(i).r[0];
       stringStream <<  ",";
       stringStream <<  population.pop.at(i).r[1];
       stringStream << ",";
       stringStream <<  population.pop.at(i).r[2];
       stringStream << ",";
       stringStream << population.pop.at(i).color;
       stringStream <<  "\n";
   }
   tempString = stringStream.str();
   buffersStreamlined.push_back(tempString);
}

void OutData::writeAndFlushBufferStreamlined() {
    std::cout << "[OutData::writeAndFlushBufferStreamlined] Writing to directory: " << OUTPUT_DIR << " (STARTED)" << std::endl;
    for (unsigned int i=0; i < buffersStreamlined.size(); i++) {
        std::ofstream outputFile;
        // outputFile.open(OUTPUT_DIR + std::to_string(i) + ".csv");
        outputFile.open(OUTPUT_DIR + std::to_string(0) + "aa.csv." + std::to_string(i), std::ofstream::trunc);
        outputFile << buffersStreamlined.at(i);
        outputFile.close();
    }
    buffersStreamlined.clear();
    std::cout << "[OutData::writeAndFlushBufferStreamlined] Writing to directory: " << OUTPUT_DIR << " (DONE)" << std::endl;
}

#endif