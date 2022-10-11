// Definitions for the functions in the PenningTrap class

#include "class_PenningTrap.hpp"
#include "Particle.hpp"

PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in){
  // Use the input variables (c0, c1) to assign values to the class memeber variables (c0_, c1_)
  //arma::vec particles = arma::vec(0);
  B0 = B0_in;
  V0 = V0_in;
  d = d_in;
}

 // void PenningTrap::add_particle(Particle p_in){
  //  particles.push_back(p_in);
 // };

  // External electric field at point r=(x,y,z)
  arma::vec PenningTrap::external_E_field(arma::vec r){ 
    double E_x = (V0*r(0))/(pow(d,2));
    double E_y = (V0*r(1))/(pow(d,2));
    double E_z = -(2*V0*r(2))/(pow(d,2));
    return arma::vec("E_x E_y E_z");
  };

  // External magnetic field at point r=(x,y,z)
  arma::vec PenningTrap::external_B_field(arma::vec r){
    double B_x = 0;
    double B_y = 0;
    double B_z = B0;
    return arma::vec("B_x B_y B_z");
  };

  // Force on particle_i from particle_j
    arma::vec PenningTrap::force_particle(int i, int j){
    double E_x_ij = k_e*particles(j).q*(particles(i).r(0)-particles(j).r(0))/pow(abs(particles(i).r(0)-particles(j).r(0)),3);
    double E_y_ij = k_e*particles(j).q*(particles(i).r(1)-particles(j).r(1))/pow(abs(particles(i).r(1)-particles(j).r(1)),3);
    double E_z_ij = k_e*particles(j).q*(particles(i).r(2)-particles(j).r(2))/pow(abs(particles(i).r(2)-particles(j).r(2)),3);
    return arma::vec("particles(i).q*E_x_ij particles(i).q*E_y_ij particles(i).q*E_z_ij");
  };  
  

    // The total force on particle_i from the external fields
    //arma::vec PenningTrap::total_force_external(int i);


    // The total force on particle_i from the other particles
    //arma::vec PenningTrap::total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles
    //arma::vec PenningTrap::total_force(int i);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    //void PenningTrap::evolve_RK4(double dt);

    // Evolve the system one time step (dt) using Forward Euler
    //void PenningTrap::evolve_forward_Euler(double dt);