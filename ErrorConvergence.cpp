
#include <armadillo>
#include <math.h>
#include "max_abs_error.hpp"


int main(){
  

  // Defining the different number of steps
  double n1 = 4000.;
  double n2 = 8000.;
  double n3 = 16000.;
  double n4 = 32000.;

  //Vector with step sizes
  arma::vec h_k = {50./n1, 50./n2, 50./n3, 50./n4};
  //Vector with number of steps
  arma::vec n_k = {n1, n2, n3, n4};

  double err_x;
  double err_y;
  double err_z;

  // Calculating error convergence rate
  for(int i=1; i<h_k.size(); i++){

    std::vector<double> abs_error = max_rel_error(n_k(i));
    std::vector<double> abs_error_prev = max_rel_error(n_k(i-1));

    err_x += 1/3*log(abs_error[0]/abs_error_prev[0])/log(h_k(i)/h_k(i-1));
    err_y += 1/3*log(abs_error[1]/abs_error_prev[1])/log(h_k(i)/h_k(i-1));
    err_z += 1/3*log(abs_error[2]/abs_error_prev[2])/log(h_k(i)/h_k(i-1));

  }

  //Printing the error convergence rates for x,y and z direction
  std::cout << err_x;
  std::cout << "#########";
  std::cout << err_y;
  std::cout << "#########";
  std::cout << err_z;


  return 1;
  }