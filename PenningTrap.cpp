// Definitions for the functions in the PenningTrap class

#include "PenningTrap.hpp"
#include "Particle.hpp"

PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in){
  // Use the input variables (c0, c1) to assign values to the class memeber variables (c0_, c1_)
  //arma::vec particles = arma::vec(0);
  B0 = B0_in;
  V0 = V0_in;
  d = d_in;
}

  void PenningTrap::add_particle(Particle p_in){
    particles.push_back(p_in);
  };

  // External electric field at point r=(x,y,z)
  arma::vec PenningTrap::external_E_field(arma::vec r){ 
    double E_x = (V0*r(0))/(pow(d,2));
    double E_y = (V0*r(1))/(pow(d,2));
    double E_z = -(2*V0*r(2))/(pow(d,2));
    arma::vec E = arma::vec(3);
    E(0) = E_x;
    E(1) = E_y;
    E(2) = E_z;
    return E;
  };

  // External magnetic field at point r=(x,y,z)
  arma::vec PenningTrap::external_B_field(arma::vec r){
    double B_x = 0;
    double B_y = 0;
    double B_z = B0;
    arma::vec B = arma::vec(3);
    B(0) = B_x;
    B(1) = B_y;
    B(2) = B_z;
    return B;
  };

  // Force on particle_i from particle_j
    arma::vec PenningTrap::force_particle(int i, int j){
    double E_x_ij = (k_e*particles[j].q*(particles[i].r(0))-particles[j].r(0))/(pow(abs(particles[i].r(0)-particles[j].r(0)),3));
    double E_y_ij = (k_e*particles[j].q*(particles[i].r(1))-particles[j].r(1))/(pow(abs(particles[i].r(1)-particles[j].r(1)),3));
    double E_z_ij = (k_e*particles[j].q*(particles[i].r(2))-particles[j].r(2))/(pow(abs(particles[i].r(2)-particles[j].r(2)),3));


    arma::vec E_force = arma::vec(3);
    E_force(0) = E_x_ij*particles[i].q;
    E_force(1) = E_y_ij*particles[i].q;
    E_force(2) = E_z_ij*particles[i].q;

    return E_force;
  };  
  

    // The total force on particle_i from the external fields
    arma::vec PenningTrap::total_force_external(int i){

    double F_x = (particles[i].q*(V0*particles[i].r(0))/(pow(d,2)))+particles[i].q*particles[i].v(1)*B0;
    double F_y = (particles[i].q*(V0*particles[i].r(1))/(pow(d,2)))-particles[i].q*particles[i].v(0)*B0; 
    double F_z = particles[i].q*(-2*V0*particles[i].r(2))/(pow(d,2));

    arma::vec F_em = arma::vec(3);
    F_em(0) = F_x;
    F_em(1) = F_y;
    F_em(2) = F_z;

    return F_em;

    }


    // The total force on particle_i from the other particles
      arma::vec PenningTrap::total_force_particles(int i){
      arma::vec E_sum = arma::vec(3).zeros();
      for(int j =0; j < particles.size(); j++){
        if(j!=i){
          E_sum += force_particle(i,j);
        }

      }
    return E_sum;
}


    // The total force on particle_i from both external fields and other particles
    arma::vec PenningTrap::total_force(int i){
      arma::vec F_tot = arma::vec(3).zeros();
      F_tot = total_force_particles(i)+total_force_external(i);

      return F_tot;
    }
    // Evolve the system one time step (dt) using Runge-Kutta 4th order

    // What abou the time dependency... (!!!)

    void PenningTrap::evolve_RK4(double dt){

      for(int i = 0; i < particles.size(); i++){
        arma::vec v_cop = particles[i].v;
        arma::vec r_cop = particles[i].r;

        arma::vec K1_r = dt*particles[i].v;
        arma::vec K1_v = (total_force(i)/particles[i].m)*dt;

        arma::vec r1 = particles[i].r + (1/2)*K1_r;
        arma::vec v1 = particles[i].v + (1/2)*K1_v;
        
        particles[i].v = v1;
        particles[i].r = r1;

        arma::vec K2_r = dt*particles[i].v;
        arma::vec K2_v = (total_force(i)/particles[i].m)*dt;

        //std::cout << K2_v;

        arma::vec r2 = particles[i].r + (1/2)*K2_r;
        arma::vec v2 = particles[i].v + (1/2)*K2_v;
        
        particles[i].v = v2;
        particles[i].r = r2;

        arma::vec K3_r = dt*particles[i].v;
        arma::vec K3_v = (total_force(i)/particles[i].m)*dt;

        //std::cout << K3_v;

        arma::vec r3 = particles[i].r + (1/2)*K3_r;
        arma::vec v3 = particles[i].v + (1/2)*K3_v;
        
        particles[i].v = v3;
        particles[i].r = r3;

        arma::vec K4_r = dt*particles[i].v;
        arma::vec K4_v = (total_force(i)/particles[i].m)*dt;
        
        std::cout << (1/6)*(K1_r + 2*K2_r + 2*K3_r + K4_r);
        std::cout << (1/6)*(K1_r + 2*K2_r + 2*K3_r + K4_r);

        arma::vec r_next = r_cop + (1/6)*(K1_r + 2*K2_r + 2*K3_r + K4_r);
        arma::vec v_next = v_cop + (1/6)*(K1_v+2*K2_v + 2*K3_v + K4_v);
        
        particles[i].v = v_next;
        particles[i].r = r_next;
      }

    }

    // Evolve the system one time step (dt) using Forward Euler
    void PenningTrap::evolve_forward_Euler(double dt){

      for(int i = 0; i < particles.size(); i++){

        arma::vec F_e = total_force(i);
        arma::vec v_next = particles[i].v + (F_e/particles[i].m)*dt;
        particles[i].v = v_next;

        arma::vec r_next = particles[i].v * dt + particles[i].r;
        particles[i].r = r_next;
      }

      //double r_x = A1*cos(-w1*t-theta)+A2*cos(-w2*t-theta)
      //double r_y = A1*sin(-w1*t-theta)+A2*sin(-w2*t-theta)


    }