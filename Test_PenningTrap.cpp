#include <armadillo>
#include "Particle.cpp"
#include "PenningTrap.cpp"
#include <math.h>

int main(){

    // Define values for the Penning trap
    double B0_in = 9.65*pow(10, 1);
    double V0_in = 9.65*pow(10, 8);
    double d_in = pow(10,4);

    // Create a Penning trap
    PenningTrap trap1 = PenningTrap(B0_in, V0_in, d_in);

    // Test that every value in the Penning trap is correct
    assert(trap1.B0=B0_in);
    assert(trap1.V0=V0_in);
    assert(trap1.d=d_in);

    // Create velocity and position for the particles
    arma::vec v1 = arma::vec("1. 2. 1.");
    arma::vec r1 = arma::vec("1. 2. 1.");

    arma::vec v2 = arma::vec("1. 2. 1.");
    arma::vec r2 = arma::vec("3. 3. 3.");

    // Create particles
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
     assert(floor(F_tot_ex(0))==202);
     assert(floor(F_tot_ex(1))==-78);
     assert(floor(F_tot_ex(2))==-20);


    // Testing function for total force on particle_i from the other particles

    arma::vec v3 = arma::vec("1. 2. 1.");
    arma::vec r3 = arma::vec("5. 5. 5.");

    Particle p3 = Particle(5., 5., r3, v3);
    trap1.add_particle(p3);

    arma::vec F_ij1 = trap1.force_particle(0,1);
    arma::vec F_ij2 = trap1.force_particle(0,2);
    arma::vec F_tot_ij_test = F_ij1 + F_ij2;

    arma::vec F_tot_ij = trap1.total_force_particles(0);
 
    assert((F_tot_ij(0))==F_tot_ij_test(0));
    assert((F_tot_ij(1))==F_tot_ij_test(1));
    assert((F_tot_ij(2))==F_tot_ij_test(2));


    // Testing the function for total force on particle

    arma::vec F_ex = trap1.total_force_external(0);
    arma::vec F_int = trap1.total_force_particles(0);
    arma::vec F_tot_test = F_ex + F_int;

    arma::vec F_tot = trap1.total_force(0);
    
    assert(F_tot(0)==F_tot_test(0));
    assert(F_tot(1)==F_tot_test(1));
    assert(F_tot(2)==F_tot_test(2));

    return 0;
}