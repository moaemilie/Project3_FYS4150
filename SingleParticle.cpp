#include <armadillo>
#include "Particle.cpp"
#include "PenningTrap.cpp"
#include <math.h>

int main(){
  
/*   // Defining values for the particels
  arma::vec r1 = arma::vec("20. 0. 20.");
  arma::vec v1 = arma::vec("0. 25. 0.");

  // Defining values for penning trap
  double B0_in = 9.65*pow(10, 1); // Rett opp
  double V0_in = 2.41*pow(10, 6); // Rett opp
  double d_in = 500.;

  double m = 40.078*1.660*pow(10,-6);
  double q = 1; 

  // Creating instsances of two particles
  Particle particle1 = Particle(q, m , r1, v1);

  // Creating penning trap
  PenningTrap trap1 = PenningTrap(B0_in, V0_in, d_in);

  // Adding particles to the penning trap
  trap1.add_particle(particle1);

  // Running the Runge kutta function 
  double n1 = 4000;
  double n2 = 8000;
  double n3 = 1600;
  double n4 = 3200;

  double TotTime = 50;
  double dt = TotTime/n1;
  int steps = n1;

  // Saving all the x-values for plotting
  std::vector<double> part1_x_nk;
  std::vector<double> part1_y_nk;
  std::vector<double> part1_z_nk;

  for(int t = 0; t < steps; t++){
      //trap1.evolve_forward_Euler(dt);  //evolve_forward_Euler;
      trap1.evolve_RK4(dt);          //evolve_RK4;
      part1_x_nk.push_back(trap1.particles[0].r(0));
      part1_y_nk.push_back(trap1.particles[0].r(1));
      part1_z_nk.push_back(trap1.particles[0].r(2));
    }

  // Write the vectors to files
  std::string filename = "SingleParticle_n1.txt";
  std::ofstream ofile;
  ofile.open(filename);
  int width = 12;
  int prec  = 4;

  // Loop over steps
  for (int i = 0; i < part1_x_n1.size(); i++)
  {
  ofile << std::setw(width) << std::setprecision(prec) << std::scientific << part1_x_n1[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part1_y_n1[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part1_z_n1[i]
          << std::endl;
  }  
  ofile.close();
 */

  // Analytical solution

  arma::vec r1 = arma::vec("20. 0. 20.");
  arma::vec v1 = arma::vec("0. 25. 0.");

  // Defining values for penning trap
  double B0 = 9.65*pow(10, 1); // Rett opp
  double V0 = 2.41*pow(10, 6); // Rett opp
  double d = 500.;
  double m = 40.078;//*1.660*pow(10,-6);
  double q = 1;


  double phi_minus = 0;
  double phi_plus = 0;
  double w0 = (q*B0)/(m);
  double wz2 = (2*q*V0)/(m*pow(d,2));
  double w_plus = (w0 + sqrt(pow(w0,2)*2*wz2))/(2);
  double w_minus = (w0 - sqrt(pow(w0,2)*2*wz2))/(2);

  std::vector<double> solution_x_analy;
  std::vector<double> solution_y_analy;
  std::vector<double> solution_z_analy;
  solution_x_analy.push_back(r1(0));
  solution_y_analy.push_back(r1(1));
  solution_z_analy.push_back(r1(2))

//Finne en måte å oppdatere v på
  for(double i = 1; i < 4000; i++){
    double t = (50./4000.)*i;
    std::cout << t;
    double A_plus = (v1(1)+(w_minus*solution_x_analy(i-1)))/(w_minus-w_plus);
    double A_minus = -(v1(1)+(w_plus*olution_x_analy(i-1)))/(w_minus-w_plus);
    double x_analy = (A_plus + cos(-w_plus*t - phi_minus)+ A_minus*cos(-w_minus*t-phi_plus));
    double y_analy = (A_minus + sin(-w_plus*t - phi_minus)+ A_minus*sin(-w_minus*t-phi_plus));
    double z_analy = r1(2)*cos(wz2*t);
    solution_x_analy.push_back(x_analy);
    solution_y_analy.push_back(y_analy);
    solution_z_analy.push_back(z_analy);
}

    // Write the vectors to files
  std::string filename = "Analytical_SingleParticle.txt";
  std::ofstream ofile;
  ofile.open(filename);
  int width = 12;
  int prec  = 4;

  // Loop over steps
  for (int i = 0; i < solution_x_analy.size(); i++)
  {
  ofile << std::setw(width) << std::setprecision(prec) << std::scientific << solution_x_analy[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << solution_y_analy[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << solution_z_analy[i]
          << std::endl; 
  }  
  ofile.close();
  
  return 1; 
}
