
#include <armadillo>
#include "Particle.cpp"
#include "PenningTrap.cpp"
#include <math.h>
#include <iostream>



int main(){

  // Set random seed
  arma::arma_rng::set_seed_random();

  // Defining values for penning trap
  double B0_in = 9.65*pow(10, 1); // Rett opp
  double V0_in = 2.41*pow(10, 6); // Rett opp
  double d_in = 500.;
  double f_in = 0.1;
  double m = 40.078;
  double q = 1; 
  double TotTime = 500.;
  double dt = 0.1;
  int steps = TotTime / dt;
  double W_steps = 0.002;  //period  

  std::vector<double> frequency;
  std::vector<double> W;

  for(double j=2; j<2.4; j=j+W_steps){

      // Creating penning trap
     PenningTrap trap_v = PenningTrap(B0_in, V0_in, d_in, f_in);

    // Adding particles to the penning trap

    for(int k = 1; k < 10; k++){

        //Filling our Penning trap with 100 random Ca+ particles
        arma::vec r = arma::vec(3).randn() * 0.1 * trap_v.d;  // random initial position
        arma::vec v = arma::vec(3).randn() * 0.1 * trap_v.d;  // random initial velocity
        Particle particle = Particle(q, m , r, v);
        trap_v.add_particle(particle);
     }

    for(int t = 1; t < steps; t++){ 
        double ev_time = t*dt;
        trap_v.evolve_RK4(dt, true, j, ev_time, true);
        //std::cout << trap_v.time_dep_V_field(j, t, true);
    } 
    double frac_left = trap_v.count_particles()/10.;
    frequency.push_back(frac_left);
    W.push_back(j);
    //std::cout << trap_v.f;
  }

  
  std::string filename2 = "Time_dep_V0_0_1_fine2.txt";
  std::ofstream ofile2;
  ofile2.open(filename2);
  int width = 12;
  int prec  = 4;
  // Loop over steps
  for (int i = 0; i < frequency.size(); i++)
  {
  ofile2 << std::setw(width) << std::setprecision(prec) << std::scientific << frequency[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << W[i]
          //<< std::setw(width) << std::setprecision(prec) << std::scientific << part2_y[i]
          //<< std::setw(width) << std::setprecision(prec) << std::scientific << part2_v_y[i]
          << std::endl;
  }  
  ofile2.close();  
  return 1;



}