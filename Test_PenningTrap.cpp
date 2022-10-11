#include <armadillo>
#include "Particle.cpp"
#include "PenningTrap.cpp"
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
    arma::vec r2 = arma::vec("3. 3. 3.");
    
    Particle p1 = Particle(1., 1., r1, v1);
    Particle p2 = Particle(1., 2., r2, v2);

    //Testing built in add_particle function

    trap1.add_particle(p1);

    assert(trap1.particles.size()==1);

    trap1.add_particle(p2);
    
    assert(trap1.particles.size()==2);

    //Testing built in External electric field function

    arma::vec r = arma::vec("2. 3. 1.");
    arma::vec E = trap1.external_E_field(r);

    assert(E(0)==19.3);
    assert(E(1)==28.95);
    assert(E(2)==-19.3);

    //Testing function for Force on particle_i from particle_j
    
    arma::vec E_force = trap1.force_particle(0,1);
    assert(floor(E_force(0))==17366);
    assert(floor(E_force(1))==277867);
    assert(floor(E_force(2))==17366);

    //Testing built in External magnetic field function

    arma::vec B = trap1.external_B_field(r);
    assert((B(0))==0);
    assert((B(1))==0);
    assert((B(2))==B0_in);

    //Testing function for total external force
    arma::vec F_tot_ex = trap1.total_force_external(0);
    assert(floor(F_tot_ex(0))==183);
    assert(floor(F_tot_ex(1))==77);
    assert(floor(F_tot_ex(2))==-20);

    // Testing function for total force on particle_i from the other particles
    
    //assert(floor(F_tot_ij(0))==17366);
    //assert(floor(F_tot_ij(1))==277867);
    //assert(floor(F_tot_ij(2))==17366);

    arma::vec v3 = arma::vec("1. 2. 1.");
    arma::vec r3 = arma::vec("5. 5. 5.");
    
    Particle p3 = Particle(5., 5., r3, v3);
    trap1.add_particle(p3);

    arma::vec F_tot_ij = trap1.total_force_particles(0); 
    assert(floor(F_tot_ij(0))==28220);
    //assert(floor(F_tot_ij(1))==329330);
    assert(floor(F_tot_ij(2))==28220);
    //std::cout << floor(F_tot_ij);

    return 0;
}