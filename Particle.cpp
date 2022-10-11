// Definitions for the functions in the Particle class


#include "Particle.hpp"

Particle::Particle(double q_in, double m_in, arma::vec r_in, arma::vec v_in)
{
  // Use the input variables (c0, c1) to assign values to the class memeber variables (c0_, c1_)
   q = q_in;
   m = m_in;
   r = r_in;
   v = v_in;
}
