// Definitions for the functions in the PenningTrap class

#include "class_PenningTrap.hpp"
#include "class_Particle.hpp"

PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
  // Use the input variables (c0, c1) to assign values to the class memeber variables (c0_, c1_)
  std::vector<Particle> particles = arma::vec(0);
  B0_in = B0_in;
  V0_in = V0_in;
  d_in = d_in;
}

  void PenningTrap::add_particle(Particle p_in){
    particles.push_back(p_in);
  };

  // External electric field at point r=(x,y,z)
  arma::vec PenningTrap::external_E_field(arma::vec r){ 
    double E_x = (V0_in*r(0))/(pow(d_in,2));
    double E_y = (V0_in*r(1))/(pow(d_in,2));
    double E_z = -(2*V0_in*r(2))/(pow(d_in,2));
    return arma::vec(E_x, E_y, E_z);
  };

  // External magnetic field at point r=(x,y,z)
  arma::vec PenningTrap::external_B_field(arma::vec r){
    double B_x = 0;
    double B_y = 0;
    double B_z = B0_in;
    return arma::vec(B_x, B_y, B_z);
  };

  // Force on particle_i from particle_j
  arma::vec PenningTrap::force_particle(int i, int j){
    double E_x_ij = k_e*particle(j).q*(particle(i).r(0)-particle(j).r(0))/pow(abs(particle(i).r(0)-particle(j).r(0)),3);
    double E_y_ij = k_e*particle(j).q*(particle(i).r(1)-particle(j).r(1))/pow(abs(particle(i).r(1)-particle(j).r(1)),3);
    double E_z_ij = k_e*particle(j).q*(particle(i).r(2)-particle(j).r(2))/pow(abs(particle(i).r(2)-particle(j).r(2)),3);
    return arma::vec(particle(i).q*E_x_ij, particle(i).q*E_y_ij, particle(i).q*E_z_ij);

        // The total force on particle_i from the external fields
    arma::vec PenningTrap::total_force_external(int i);


    // The total force on particle_i from the other particles
    arma::vec PenningTrap::total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles
    arma::vec PenningTrap::total_force(int i);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void PenningTrap::evolve_RK4(double dt);

    // Evolve the system one time step (dt) using Forward Euler
    void PenningTrap::evolve_forward_Euler(double dt);

  }