#ifndef AAVECTOR_H
#define AAVECTOR_H

/* This is where vector math lives */
#include <vector>
#include <iostream>
#include <cmath>

#include <cfenv>


namespace aa {

// vector dot product
inline double dot(std::vector<double> a, std::vector<double> b) {
    if ((a.size() == 2 || a.size() == 3) && a.size() == b.size()) {
        double sum = 0;
        for (unsigned int i=0; i < a.size(); i++) {
            sum += a[i]*b[i];
        }
        return sum;
    } else {
        std::cout << "<<ERROR>> [aa::dot] vector size error" << std::endl;
        return 0.0;
    }     
}

// vector componentwise product
inline std::vector<double> compprod(std::vector<double> a, std::vector<double> b) {
     if (a.size() == 2 || a.size() == 3) {
        std::vector<double> ret;
        for (unsigned int i=0; i < a.size(); i++) {
            ret.push_back(a[i] * b[i]);
        }
        return ret;
    } else {
        std::cout << "<<ERROR>> [aa::dot] vector size error" << std::endl;
        return {0.0, 0.0, 0.0};
    }     
}

// multiply a vector by a constant
inline std::vector<double> cmul(double constant, std::vector<double> a) {
    if (a.size() == 2 || a.size() == 3) {
        std::vector<double> ret;
        for (unsigned int i=0; i < a.size(); i++) {
            ret.push_back(constant * a[i]);
        }
        return ret;
    } else {
        std::cout << "<<ERROR>> [aa::dot] vector size error" << std::endl;
        return {0.0, 0.0, 0.0};
    }     
}

// subtract a constant from a vector
inline std::vector<double> csub(double constant, std::vector<double> a) {
    if (a.size() == 2 || a.size() == 3) {
        std::vector<double> ret;
        for (unsigned int i=0; i < a.size(); i++) {
            ret.push_back(a[i] - constant);
        }
        return ret;
    } else {
        std::cout << "<<ERROR>> [aa::dot] vector size error" << std::endl;
        return {0.0, 0.0, 0.0};
    }     
}

// add a constant to a vector
inline std::vector<double> cadd(double constant, std::vector<double> a) {
    if (a.size() == 2 || a.size() == 3) {
        std::vector<double> ret;
        for (unsigned int i=0; i < a.size(); i++) {
            ret.push_back(a[i] + constant);
        }
        return ret;
    } else {
        std::cout << "<<ERROR>> [aa::dot] vector size error" << std::endl;
        return {0.0, 0.0, 0.0};
    }     
}

// subtract vector b from vector a
inline std::vector<double> sub(std::vector<double> a, std::vector<double> b) {
    std::vector<double> ret;
    if (a.size() == b.size()) {
        for (unsigned int i=0; i < a.size(); i++) {
            ret.push_back(a[i] - b[i]);
        }
        return ret;
    }
    else {
        std::cout << "<<ERROR>> [aa::sub] vector size error" << std::endl;
        return {0.0, 0.0, 0.0};
    }
}

// add two vectors
inline std::vector<double> add(std::vector<double> a, std::vector<double> b) {
    std::vector<double> ret;
    if (a.size() == b.size()) {
        for (unsigned int i=0; i < a.size(); i++) {
            ret.push_back(a[i] + b[i]);
        }
        return ret;
    }
    else {
        std::cout << "<<ERROR>> [aa::sub] vector size error" << std::endl;
        return {0.0, 0.0, 0.0};
    }
}

// magnitude of a vector
inline double mag(std::vector<double> a) {
    if (a.size() == 0) {
        return 0;
    } else {
        double sum = 0.0;
        for (auto i : a) {
            sum += i*i;
        }
        return sqrt(sum);
    }
}

// magnitude of a vector squared
inline double mag2(std::vector<double> a) {
    if (a.size() == 0) {
        return 0;
    } else {
        double sum = 0.0;
        for (auto i : a) {
            sum += i*i;
        }
        return sum;
    }
}

};

#endif