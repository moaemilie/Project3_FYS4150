#ifndef __filename_hpp__
#define __filename_hpp__

#include <armadillo>
#include <math.h>
#include "Particle.hpp"
#include <iostream>

class PenningTrap{

    private:

    public:
        double B0;
        double V0;
        double d;
        double f;
        std::vector<Particle> particles;
        //arma::vec particles;
        double k_e = 1.38935333*pow(10,5);

        PenningTrap(){};
        // Constructor
        PenningTrap(double B0_in, double V0_in, double d_in, double f_in);

        // Add a particle to the trap
        void add_particle(Particle p_in);
        
        // Time depended V-field where the time-dependency can be turned on and off
        double time_dep_V_field(double w_v, double t, bool v_field);

        // External electric field at point r=(x,y,z)
        arma::vec external_E_field(arma::vec r, double w_v, double t, bool v_field);

        // External magnetic field at point r=(x,y,z)
        arma::vec external_B_field(arma::vec r);

        // Force on particle_i from particle_j
        arma::vec force_particle(int i, int j);

        // The total force on particle_i from the external fields
        arma::vec total_force_external(int i, double w_v, double t, bool v_field);

        // The total force on particle_i from the other particles
        arma::vec total_force_particles(int i);

        // The total force on particle_i from both external fields and other particles
        arma::vec total_force(int i, bool Coulomb, double w_v, double t, bool v_field);

        // Evolve the system one time step (dt) using Runge-Kutta 4th order
        void evolve_RK4(double dt, bool Coulomb, double w_v, double t, bool v_field);

        // Evolve the system one time step (dt) using Forward Euler
        void evolve_forward_Euler(double dt, bool Coulomb, double w_v, double t, bool v_field);

        //Counts the number of particles that are left in the trap
        double count_particles();

};


#endif