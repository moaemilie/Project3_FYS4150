#include <armadillo>
#include "Particle.cpp"
#include "Copy_Penning_Trap.cpp"
#include <math.h>
//#include "max_abs_error.hpp"

int main(){
  
  // Defining values for the particel
  arma::vec r1 = arma::vec("20. 0. 20.");
  arma::vec v1 = arma::vec("0. 25. 0.");
  double m = 40.078;
  double q = 1; 

  // Defining values for penning trap
  double B0_in = 9.65*pow(10, 1); 
  double V0_in = 2.41*pow(10, 6); 
  double d_in = 500.;
  double w_v = 2.5;

  // Creating instsances of a particel
  Particle particle1 = Particle(q, m , r1, v1);

  // Creating a penning trap
  PenningTrap trap1 = PenningTrap(B0_in, V0_in, d_in);

  // Adding particles to the penning trap
  trap1.add_particle(particle1);

  // Defining values for the different dt-s
  double n1 = 4000;
  double n2 = 8000;
  double n3 = 16000;
  double n4 = 32000;

  double TotTime = 50;
  double dt = TotTime/n4;

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

  // Make a place to save the analytical solutions
  std::vector<double> solution_x_analy;
  std::vector<double> solution_y_analy;
  std::vector<double> solution_z_analy;

  // Make a place to save the relative error
  std::vector<double> Error_rel_x_nk;
  std::vector<double> Error_rel_y_nk;
  std::vector<double> Error_rel_z_nk;

    for(int i = 1; i < n4; i++){
      //trap1.evolve_RK4(dt, false);          //evolve with Runge Kutta, this was done before the time modification of the PenningTrap class
      trap1.evolve_forward_Euler(dt, false); //evolve with Euler, this was done before the time modification of the PenningTrap class
      part1_x_nk.push_back(trap1.particles[0].r(0));
      part1_y_nk.push_back(trap1.particles[0].r(1));
      part1_z_nk.push_back(trap1.particles[0].r(2));

      // Calculate the analytical solution
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
      Error_rel_x_nk.push_back((trap1.particles[0].r(0) - x_analy)/(trap1.particles[0].r(0)));
      Error_rel_y_nk.push_back((trap1.particles[0].r(1) - y_analy)/(trap1.particles[0].r(1)));
      Error_rel_z_nk.push_back((trap1.particles[0].r(2) - z_analy)/(trap1.particles[0].r(2)));
    }


  // Write the vectors to files
  std::string filename = "Rel_error_n4_COPY_EUL.txt";
  std::ofstream ofile;
  ofile.open(filename);
  int width = 12;
  int prec  = 4;

  // Loop over steps
  for (int i = 0; i < Error_rel_x_nk.size(); i++){
  ofile << std::setw(width) << std::setprecision(prec) << std::scientific << Error_rel_x_nk[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << Error_rel_y_nk[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << Error_rel_z_nk[i]
          << std::endl; 
  }  
  ofile.close();   

return 1;
}
