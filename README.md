# Project3_FYS4150
# Project3_FYS4150

## Our files

- In the source file Particle.cpp we have our particle class. In this file we have also included the header file Particel.hpp containg the class declaration
  for the particle class

- In the source file PenningTrap.cpp we have our Penning Trap class. In this class we have included the particle class from Particle.cpp and 
  the header file PenningTrap.hpp containing the class declaration for our Penning trap class 
  
- In the sourcefile Task8.cpp we use the Penning Trap class from PenningTrap.cpp and evolve the penning trap in time by using nummeric functions imbedded 
  in the Penning trap class.
  
- In RelativeError.cpp we have included the Penning trap class and here we run the nummeric functions from the penning trap class and find the relative 
  error between the nummeric and analytical solutions.
  
- In max_abs_error.hpp we have made a function that first finds the absolute error between the nummeric and the analytical solution and the finds the max
  absolute error
  
- In the sourcefile ErrorConvergence.cpp we include the the header file max_abs_error.cpp and use the function to calculate the error convergence.

- In Task9.cpp we include the PenningTrap.cpp file and simulate the system with a time dependent electric field and calculate the fraction of particles left in the trap after a certain amount of time.

- Our tests
    - Test_PenningTrap.cpp
    - Test_Particle.cpp

## How to Compile and run the files

- Task8.cpp, RealtiveError.cpp and ErrorConvergence.cpp were made and run before the time dependent electric potential was added to the Penning trap class. To run these files with the new version, the if-loop in the external E-field and B-field has to be commented. Also the extra inputs has to be added to the Euler and Runge Kutta function.
- To run Task8.cpp the parts of the file you want to run has to be uncommented. The parts should be run both with and without coloumb force. This can be turned on by putting true as second argument in the Runge Kutta or Euler function. Each part should be run for both the Euler and the Runge Kutta function. Each part saves a file with its data.
- RelativeError.cpp should be run for the different step sizes, n1, n2, n3 and n4, and the names of the files that the data is written to has to be changed for each stepsize.
- Error Convergence should be compiled and run to produce the error convergence rates
- Task9 has to be compiled and run for the frequencies: f=0.1, f=0.4, f=0.7.  The simulation should also be run for a range of frequencies between 2 and 2.4 and with a step size of 0.002. The simulation should also be run for the different amplitudes in this range. Then the name of the files that the data is written to should be changed for each run
- To test the Penning trap class run and compile Test_PenningTrap.cpp
- To test the Particle class run and compile the Test_Particle.cpp
