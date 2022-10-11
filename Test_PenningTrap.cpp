#include <armadillo>
#include "Particle.cpp"
#include "class_PenningTrap.cpp"
#include <math.h>

int main(){

    double B0_in = 9.65*pow(10, 1);
    double V0_in = 9.65*pow(10, 8);
    double d_in = pow(10,4);

    PenningTrap trap1 = PenningTrap(B0_in, V0_in, d_in);

    assert(trap1.B0=B0_in);
    assert(trap1.V0=V0_in);
    assert(trap1.d=d_in);


    arma::vec v1 = arma::vec("1. 2. 1.");
    arma::vec r1 = arma::vec("1. 2. 1.");

    arma::vec v2 = arma::vec("1. 2. 1.");
    arma::vec r2 = arma::vec("1. 2. 1.");
    
    Particle p1 = Particle(1., 1., r1, v1);
    Particle p2 = Particle(1., 2., r2, v2);

    return 0;
}