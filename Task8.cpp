#include <armadillo>
#include "Particle.cpp"
#include "PenningTrap.cpp"
#include <math.h>

    arma::vec v1 = arma::vec("0. 25. 0.");
    arma::vec r1 = arma::vec("20.*pow(10,-6) 0. 20.*pow(10,-6)");

    arma::vec v2 = arma::vec("0. 40. 5.");
    arma::vec r2 = arma::vec("25.*pow(10,-6) 25.*pow(10,-6) 0.");
    
    Particle particle1 = Particle(1., 1., r1, v1);
    Particle particle2 = Particle(1., 1., r2, v2);

    int TotTime = 50*pow(10,-6);
    for(int t = 0; t < TotTime){

      }
