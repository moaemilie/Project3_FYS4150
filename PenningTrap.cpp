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

    double E_x = (V0*r(0))/(pow(d,2));  // Her virker det som at den runder av når tallene blir små... 0.00016 blir 0.0002 (ikke bra)
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

      double E_x_ij = (k_e*particles[j].q*(((particles[i].r(0))-(particles[j].r(0)))/(pow((abs(particles[i].r(0)-particles[j].r(0))),3))));      
      double E_y_ij = (k_e*particles[j].q*(((particles[i].r(1))-(particles[j].r(1))))/(pow((abs(particles[i].r(1)-particles[j].r(1))),3)));
      double E_z_ij = (k_e*particles[j].q*(((particles[i].r(2))-(particles[j].r(2))))/(pow((abs(particles[i].r(2)-particles[j].r(2))),3)));

      arma::vec E_force = arma::vec(3);
      E_force(0) = E_x_ij * particles[i].q;
      E_force(1) = E_y_ij * particles[i].q;
      E_force(2) = E_z_ij * particles[i].q;

      return E_force;
    };  
  

    // The total force on particle_i from the external fields

    // Skulle kanskje bruke våre tidligere funkjsoner her? ... (!!!)

    arma::vec PenningTrap::total_force_external(int i){
      arma::vec E_field = external_E_field(particles[i].r);

      double F_x = (particles[i].q*E_field(0))+particles[i].q*particles[i].v(1)*B0;
      double F_y = (particles[i].q*E_field(1))-particles[i].q*particles[i].v(0)*B0; 
      double F_z = particles[i].q*E_field(2);

      arma::vec F_em = arma::vec(3);
      F_em(0) = F_x;
      F_em(1) = F_y;
      F_em(2) = F_z;

      return F_em;

    }


    // The total force on particle_i from the other particles
    arma::vec PenningTrap::total_force_particles(int i){
      arma::vec E_sum = arma::vec(3).zeros();
      for(int j = 0; j < particles.size(); j++){
        if(j!=i){
          //E_sum += force_particle(i,j);
          E_sum(0) += force_particle(i,j)(0);
          E_sum(1) += force_particle(i,j)(1);
          E_sum(2) += force_particle(i,j)(2);
        }
      }
      return E_sum;
  }


    // The total force on particle_i from both external fields and other particles
  arma::vec PenningTrap::total_force(int i){
    arma::vec F_tot = arma::vec(3).zeros();
    //F_tot = total_force_particles(i)+total_force_external(i); // Force with interaction
    F_tot = total_force_external(i); // Force without interaction
    return F_tot;
  }

    // Evolve the system one time step (dt) using Runge-Kutta 4th order

    // What abou the time dependency... (!!!)

  void PenningTrap::evolve_RK4(double dt){

    std::vector<Particle> init_part;
    std::vector<arma::mat> K1_v;
    std::vector<arma::mat> K1_r;
    std::vector<arma::mat> K2_v;
    std::vector<arma::mat> K2_r;
    std::vector<arma::mat> K3_v;
    std::vector<arma::mat> K3_r;
    std::vector<arma::mat> K4_v;
    std::vector<arma::mat> K4_r;
    
    // Make a copy of all the particle velocities and positions
    for(int i = 0; i < particles.size(); i++){

      init_part.push_back(particles[i]);

      K1_r.push_back(dt*init_part[i].v);
      K1_v.push_back((total_force(i)/init_part[i].m)*dt);

      particles[i].r = init_part[i].r + (1./2.)*K1_r[i]; 
      particles[i].v = init_part[i].v + (1./2.)*K1_v[i];
      }

    for(int i = 0; i < particles.size(); i++){

      K2_r.push_back(dt*particles[i].v);
      K2_v.push_back((total_force(i)/particles[i].m)*dt);

      // Calculate r2 and v2 and updating the particle velocity and position
      particles[i].v = init_part[i].v + (1./2.)*K2_v[i];
      particles[i].r = init_part[i].r + (1./2.)*K2_r[i]; 

    } 

    for(int i = 0; i < particles.size(); i++){
      K3_r.push_back(dt*particles[i].v);
      K3_v.push_back((total_force(i)/particles[i].m)*dt);

      // Calculate r3 and v3 and updating the particle velocity and position
      particles[i].v = init_part[i].v + K3_v[i];
      particles[i].r = init_part[i].r + K3_r[i];
    }

    for(int i = 0; i < particles.size(); i++){
      K4_r.push_back(dt*particles[i].v);
      K4_v.push_back((total_force(i)/particles[i].m)*dt);
      
      // Uupdating the particle velocity and position
      particles[i].v = init_part[i].v + (1./6.)*(K1_v[i] + 2*K2_v[i] + 2*K3_v[i] + K4_v[i]);
      particles[i].r = init_part[i].r + (1./6.)*(K1_r[i] + 2*K2_r[i] + 2*K3_r[i] + K4_r[i]);
    }
  }

    // Evolve the system one time step (dt) using Forward Euler
  void PenningTrap::evolve_forward_Euler(double dt){

    for(int i = 0; i < particles.size(); i++){
      arma::vec v_next = particles[i].v + ((total_force(i))/(particles[i].m))*dt;
      arma::vec r_next = particles[i].r + (particles[i].v*dt);

      particles[i].v = v_next;
      particles[i].r = r_next;
    } 
  }