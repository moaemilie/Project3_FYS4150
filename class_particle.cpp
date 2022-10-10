// Definitions for the functions in the Particle class


#include "class_particle.hpp"

Particle::Particle(bool q, bool m, arma::vec r, arma::vec v)
{
  // Use the input variables (c0, c1) to assign values to the class memeber variables (c0_, c1_)
   q = q;
   m = m;
   r = r;
   v = v;
}