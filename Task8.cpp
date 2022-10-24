#include <armadillo>
#include "Particle.cpp"
#include "PenningTrap.cpp"
#include <math.h>

int main(){
  
  // Defining values for the particels
  arma::vec r1 = arma::vec("20. 0. 20.");
  arma::vec v1 = arma::vec("0. 25. 0.");

  // Defining values for penning trap
  double B0_in = 9.65*pow(10, 1); // Rett opp
  double V0_in = 2.41*pow(10, 6); // Rett opp
  double d_in = 500.;

  double m = 40.078;
  double q = 1; 

  // Creating instsances of two particles
  Particle particle1 = Particle(q, m , r1, v1);

  // Creating penning trap
  PenningTrap trap1 = PenningTrap(B0_in, V0_in, d_in);

  // Adding particles to the penning trap
  trap1.add_particle(particle1);


 
  // Running the Runge kutta function 
  double TotTime = 50.;
  double dt = 0.1;
  int steps = TotTime / dt;

  // Saving all the values for plotting
  /*std::vector<double> part1_x;
  std::vector<double> part1_y;
  std::vector<double> part1_z;

  // Add initial values
  part1_x.push_back(trap1.particles[0].r(0));
  part1_y.push_back(trap1.particles[0].r(1));
  part1_z.push_back(trap1.particles[0].r(2));

   for(int t = 1; t < steps; t++){
      trap1.evolve_RK4(dt);          //evolve_RK4, this was done before the time modification of the PenningTrap class
      part1_x.push_back(trap1.particles[0].r(0));
      part1_y.push_back(trap1.particles[0].r(1));
      part1_z.push_back(trap1.particles[0].r(2));
    }


// Write the vectors to files
  std::string filename = "Particle_z_RK4.txt";
  std::ofstream ofile;
  ofile.open(filename);
  int width = 12;
  int prec  = 4;
  // Loop over steps
  for (int i = 0; i < part1_x.size(); i++)
  {
  ofile << std::setw(width) << std::setprecision(prec) << std::scientific << part1_x[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part1_y[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part1_z[i]
          << std::endl;
  }  
  ofile.close(); */

//#####################################################################################

// Adding another particle

  arma::vec r2 = arma::vec("25. 25. 0.");
  arma::vec v2 = arma::vec("0. 40. 5.");
  Particle particle2 = Particle(q, m, r2, v2);

  trap1.add_particle(particle2);

//#####################################################################################

// Getting the x and y values for the two particles with and withount interaction
 /* 
  // Saving all the values for plotting
  std::vector<double> part1_x;
  std::vector<double> part1_y;
  std::vector<double> part2_x;
  std::vector<double> part2_y;

// Add initial values
  part1_x.push_back(trap1.particles[0].r(0));
  part1_y.push_back(trap1.particles[0].r(1));
  part2_x.push_back(trap1.particles[1].r(0));
  part2_y.push_back(trap1.particles[1].r(1));

  // Running the Runge kutta function 
  for(int t = 1; t < steps; t++){
      trap1.evolve_RK4(dt, false);          //evolve_RK4, this was done before the time modification of the PenningTrap class
      part1_x.push_back(trap1.particles[0].r(0));
      part1_y.push_back(trap1.particles[0].r(1));
      part2_x.push_back(trap1.particles[1].r(0));
      part2_y.push_back(trap1.particles[0].r(1));
    } 
  
  // Write the vectors to files
  std::string filename = "2Particle_RK4_NO_INTER.txt";
  std::ofstream ofile;
  ofile.open(filename);
  int width = 12;
  int prec  = 4;
  // Loop over steps
  for (int i = 0; i < part1_x.size(); i++)
  {
  ofile << std::setw(width) << std::setprecision(prec) << std::scientific << part1_x[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part1_y[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part2_x[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part2_y[i]
          << std::endl;
  }  
  ofile.close();  */

  //#####################################################################################
  // Getting the x, v_x, y and v_y values for the two particles with and withount interaction
/* 
  std::vector<double> part1_x;
  std::vector<double> part1_v_x;
  std::vector<double> part1_y;
  std::vector<double> part1_v_y;
  std::vector<double> part2_x;
  std::vector<double> part2_v_x;
  std::vector<double> part2_y;
  std::vector<double> part2_v_y;

  // Add initial values
  part1_x.push_back(trap1.particles[0].r(0));
  part1_v_x.push_back(trap1.particles[0].v(0));
  part1_y.push_back(trap1.particles[0].r(1));
  part1_v_y.push_back(trap1.particles[0].v(1));

  part2_x.push_back(trap1.particles[1].r(0));
  part2_v_x.push_back(trap1.particles[1].v(0));
  part2_y.push_back(trap1.particles[1].r(1));
  part2_v_y.push_back(trap1.particles[1].v(1));

  // Running the Runge kutta function 
  for(int t = 1; t < steps; t++){
      trap1.evolve_RK4(dt, false);          //evolve_RK4, this was done before the time modification of the PenningTrap class

      part1_x.push_back(trap1.particles[0].r(0));
      part1_v_x.push_back(trap1.particles[0].v(0));
      part1_y.push_back(trap1.particles[0].r(1));
      part1_v_y.push_back(trap1.particles[0].v(1));

      part2_x.push_back(trap1.particles[1].r(0));
      part2_v_x.push_back(trap1.particles[1].v(0));
      part2_y.push_back(trap1.particles[1].r(1));
      part2_v_y.push_back(trap1.particles[1].v(1));
    } 
  
  // Write the vectors to files
  std::string filename = "Particle1_RK4_NO_INTER_speed.txt";
  std::ofstream ofile;
  ofile.open(filename);
  int width = 12;
  int prec  = 4;
  // Loop over steps
  for (int i = 0; i < part1_x.size(); i++)
  {
  ofile << std::setw(width) << std::setprecision(prec) << std::scientific << part1_x[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part1_v_x[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part1_y[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part1_v_y[i]
          << std::endl;
  }  
  ofile.close();

    // Write the vectors to files
  std::string filename2 = "Particle2_RK4_NO_INTER_speed.txt";
  std::ofstream ofile2;
  ofile2.open(filename2);
  // Loop over steps
  for (int i = 0; i < part2_x.size(); i++)
  {
  ofile2 << std::setw(width) << std::setprecision(prec) << std::scientific << part2_x[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part2_v_x[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part2_y[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part2_v_y[i]
          << std::endl;
  }  
  ofile2.close(); */  
 
  //#####################################################################################
  // Getting the x, y, z values for the two particles with and withount interaction
 
  std::vector<double> part1_x;
  std::vector<double> part1_y;
  std::vector<double> part1_z;
  std::vector<double> part2_x;
  std::vector<double> part2_y;
  std::vector<double> part2_z;


  // Add initial values
  part1_x.push_back(trap1.particles[0].r(0));
  part1_y.push_back(trap1.particles[0].v(0));
  part1_z.push_back(trap1.particles[0].r(1));

  part2_x.push_back(trap1.particles[1].r(0));
  part2_y.push_back(trap1.particles[1].v(0));
  part2_z.push_back(trap1.particles[1].r(1));


  // Running the Runge kutta function 
  for(int t = 1; t < steps; t++){
      trap1.evolve_RK4(dt, false);          //evolve_RK4, this was done before the time modification of the PenningTrap class

      part1_x.push_back(trap1.particles[0].r(0));
      part1_y.push_back(trap1.particles[0].v(0));
      part1_z.push_back(trap1.particles[0].r(1));

      part2_x.push_back(trap1.particles[1].r(0));
      part2_y.push_back(trap1.particles[1].v(0));
      part2_z.push_back(trap1.particles[1].r(1));
    } 
  
  // Write the vectors to files
  std::string filename = "Particles_NO_INTER.txt";
  std::ofstream ofile;
  ofile.open(filename);
  int width = 12;
  int prec  = 4;
  // Loop over steps
  for (int i = 0; i < part1_x.size(); i++)
  {
  ofile << std::setw(width) << std::setprecision(prec) << std::scientific << part1_x[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part1_y[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part1_z[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part2_x[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part2_y[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << part2_z[i]
          << std::endl;
  }  
  ofile.close(); 

  return 1; 
}
