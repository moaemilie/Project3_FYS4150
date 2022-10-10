#ifndef __filename2_hpp__
#define __filename2_hpp__

#include <armadillo>
#include "class_Particle.hpp"

class Particle{
private:

public:
    bool q;
    bool m;
    arma::vec r;
    arma::vec v;

    Particle(){};

    Particle(bool q, bool m, arma::vec r, arma::vec v);

};

#endif