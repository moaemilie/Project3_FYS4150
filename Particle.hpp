#ifndef __filename2_hpp__
#define __filename2_hpp__

#include <armadillo>

class Particle
{
    private:

    public:

        Particle(double q_in, double m_in, arma::vec r_in, arma::vec v_in);

        double q;
        double m;
        arma::vec r;
        arma::vec v;
};

#endif