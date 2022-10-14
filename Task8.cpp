#include <armadillo>
#include "Particle.cpp"
#include "PenningTrap.cpp"
#include <math.h>

int main(){
  
  // Defining values for the particels
  arma::vec r1 = arma::vec("20. 0. 20.");
  arma::vec v1 = arma::vec("0. 25. 0.");

  arma::vec r2 = arma::vec("25. 25. 0.");
  arma::vec v2 = arma::vec("0. 40. 5.");

  // Defining values for penning trap
  double B0_in = 9.65*pow(10, 1); // Rett opp
  double V0_in = 2.41*pow(10, 6); // Rett opp
  double d_in = 500.;

  double m = 40.078*1.660*pow(10,-6);
  //double q = 1*pow(10,6); //1.602*pow(10,-19);
  double q = 1; //(V0_in*4*m)/(pow(B0_in,2)*pow(d_in,2));
  std::cout << q;

  // Creating instsances of two particles
  Particle particle1 = Particle(q, m , r1, v1);
  Particle particle2 = Particle(q, m, r2, v2);

  // Creating penning trap
  PenningTrap trap1 = PenningTrap(B0_in, V0_in, d_in);

  // Adding particles to the penning trap
  trap1.add_particle(particle1);
  trap1.add_particle(particle2);


  //double dt1 = 0.000001;
  //trap1.evolve_forward_Euler(dt1);


  // Running the Runge kutta function 
  double TotTime = 0.000050;//*pow(10,-6);
  double dt = 0.000001;
  int steps = TotTime / dt;

  // Saving all the z-values for plotting
  std::vector<double> part1_x;
  std::vector<double> part2_x;


  for(int t = 0; t < steps; t++){
      //trap1.evolve_forward_Euler(dt);  //evolve_forward_Euler;
      trap1.evolve_RK4(dt);          //evolve_RK4;
      part1_x.push_back(trap1.particles[0].r(2));
      part2_x.push_back(trap1.particles[1].r(2));
    }



// Write the vectors to files
  std::string filename = "Particle_x_RK4.txt";
  std::ofstream ofile;
  ofile.open(filename);
  int width = 12;
  int prec  = 4;

  // Loop over steps
  for (int i = 0; i < part1_x.size(); i++)
  {
  ofile << std::setw(width) << std::setprecision(prec) << std::scientific << part1_x[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part2_x[i]
          << std::endl;
  }  
  ofile.close();

  return 1; 
}
