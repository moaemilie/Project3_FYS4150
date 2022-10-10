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

  void add_particle(Particle p_in){
    particles.push_back(p_in);
  };

  // External electric field at point r=(x,y,z)
  arma::vec external_E_field(arma::vec r){ 
    bool E_x = (V0_in*r(0))/(pow(d_in,2));
    bool E_y = (V0_in*r(1))/(pow(d_in,2));
    bool E_z = -(2*V0_in*r(2))/(pow(d_in,2));
    return arma::vec(E_x, E_y, E_z);
  };

  // External magnetic field at point r=(x,y,z)
  arma::vec external_B_field(arma::vec r){
    bool B_x = 0;
    bool B_y = 0;
    bool B_z = B0_in;
    return arma::vec(B_x, B_y, B_z);
  };

  // Force on particle_i from particle_j
  arma::vec force_particle(int i, int j){
    bool E_x_ij = ;
    bool E_y_ij = ;
    bool E_z_ij = ;
    return arma::vec(E_x, E_y, E_z);
  }