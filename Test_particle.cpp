
#include <armadillo>
#include "Particle.hpp"
#include "Particle.cpp"

int main(){

    arma::vec v1 = arma::vec("1. 2. 1.");
    arma::vec r1 = arma::vec("1. 2. 1.");
    double n = 1.0;
    double h = 1.0;

    arma::vec v2 = arma::vec("1. 2. 1.");
    arma::vec r2 = arma::vec("1. 2. 1.");
    
    Particle p1 = Particle(n, h, r1, v1);
    Particle p2 = Particle(1., 2., r2, v2);

    assert(p1.v(0)=1.);
    assert(p1.r(0)=1.);
    assert(p2.v(0)=1.);
    assert(p2.r(0)=1.);

    return 0;
}