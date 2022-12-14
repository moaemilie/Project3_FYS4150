#include <armadillo>
#include "Particle.cpp"
#include "Copy_Penning_Trap.cpp"
#include <math.h>

std::vector<double> max_rel_error(double n_k);

std::vector<double> max_rel_error(double n_k){
  
  // Defining values for the particel
  arma::vec r1 = arma::vec("20. 0. 20.");
  arma::vec v1 = arma::vec("0. 25. 0.");

  // Defining values for penning trap
  double B0_in = 9.65*pow(10, 1); 
  double V0_in = 2.41*pow(10, 6); 
  double d_in = 500.;
  double w_v = 2.5;
  double m = 40.078;
  double q = 1.; 

  // Creating instsances of a particel
  Particle particle1 = Particle(q, m , r1, v1);

  // Creating a penning trap
  PenningTrap trap1 = PenningTrap(B0_in, V0_in, d_in);

  // Adding particles to the penning trap
  trap1.add_particle(particle1);

  // Defining values for the different dt-s

  double TotTime = 50.;
  double dt = TotTime/n_k;

  // Saving all the x,y,z-values for plotting
  std::vector<double> part1_x_nk;
  std::vector<double> part1_y_nk;
  std::vector<double> part1_z_nk;


  // Analytical solution

  arma::vec r1_analy = arma::vec("20. 0. 20.");
  arma::vec v1_analy = arma::vec("0. 25. 0.");

  double phi_minus = 0;
  double phi_plus = 0;
  double w0 = (q*B0_in)/(m);
  double wz2 = (2.*q*V0_in)/(m*pow(d_in,2));
  double w_plus = (w0 + (sqrt(pow(w0,2)-2.*wz2)))/(2.);
  double w_minus = (w0 - (sqrt(pow(w0,2)-2.*wz2)))/(2.);

  std::vector<double> solution_x_analy;
  std::vector<double> solution_y_analy;
  std::vector<double> solution_z_analy;

  // Define the relative error
  std::vector<double> Error_rel_x_nk;
  std::vector<double> Error_rel_y_nk;
  std::vector<double> Error_rel_z_nk;

    for(double i = 1; i < n_k; i++){
      trap1.evolve_forward_Euler(dt, true);          //evolve_RK4;
      part1_x_nk.push_back(trap1.particles[0].r(0));
      part1_y_nk.push_back(trap1.particles[0].r(1));
      part1_z_nk.push_back(trap1.particles[0].r(2));

      double t = dt*i;
      double A_plus = (v1_analy(1)+(w_minus*r1(0)))/(w_minus-w_plus);
      double A_minus = (-v1_analy(1)-(w_plus*r1(0)))/(w_minus-w_plus);
    
      double y_analy = (-A_plus*sin(w_plus*t + phi_plus)) - (A_minus*sin(w_minus*t + phi_minus));
      double x_analy = (A_plus*cos(w_plus*t + phi_plus)) + (A_minus*cos(w_minus*t + phi_minus));
      double z_analy = (r1_analy(2))*cos(sqrt(wz2)*t);

      solution_x_analy.push_back(x_analy);
      solution_y_analy.push_back(y_analy);
      solution_z_analy.push_back(z_analy);

      // Calculate the relative error
      Error_rel_x_nk.push_back(trap1.particles[0].r(0) - x_analy);
      Error_rel_y_nk.push_back(trap1.particles[0].r(1) - y_analy);
      Error_rel_z_nk.push_back(trap1.particles[0].r(2) - z_analy);
    }


  // Calculating the max relative error in x-direction
  std::vector<double> max_rel_values;
  double max_valx = Error_rel_x_nk[0];

  for(int i = 1; i < Error_rel_x_nk.size(); i++){

    if(max_valx < Error_rel_x_nk[i]){
      
      max_valx=Error_rel_x_nk[i];

    }
  }
  max_rel_values.push_back(max_valx);

 // Calculating the max relative error in y-direction
  double max_valy = Error_rel_y_nk[0];

  for(int i = 1; i < Error_rel_y_nk.size(); i++){

    if(max_valx < Error_rel_y_nk[i]){
      
      max_valx = Error_rel_y_nk[i];

    }
  }
  max_rel_values.push_back(max_valy);

 // Calculating the max relative error in z-direction
  double max_valz = Error_rel_z_nk[0];

  for(int i = 1; i < Error_rel_z_nk.size(); i++){

    if(max_valx < Error_rel_z_nk[i]){
      
      max_valx = Error_rel_z_nk[i];

    }
  }
  max_rel_values.push_back(max_valz);
  

  
  return max_rel_values; 
}
