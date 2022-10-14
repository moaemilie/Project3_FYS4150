#include <armadillo>
#include "Particle.cpp"
#include "PenningTrap.cpp"
#include <math.h>

int main(){
    arma::vec r1 = arma::vec("2. 0. 2.");
    arma::vec v1 = arma::vec("0. 1. 0.");

    arma::vec r2 = arma::vec("3. 3. 0.");
    arma::vec v2 = arma::vec("0. 4. 5.");

    // Defining values for penning trap
    double B0_in = 10; // Rett opp
    double V0_in = 20; // Rett opp
    double d_in = 500.;

    double m = 10;
    double q = 1;

    Particle particle1 = Particle(q, m , r1, v1);
    Particle particle2 = Particle(q, m, r2, v2);

    PenningTrap trap1 = PenningTrap(B0_in, V0_in, d_in);
    trap1.add_particle(particle1);
    trap1.add_particle(particle2);

    arma::vec E_ex_1 = trap1.total_force(0);
    arma::vec E_ex_2 = trap1.force_particle(0,1);
    std::cout << E_ex_1;
    std::cout << "#";
    //std::cout << E_ex_2;
}