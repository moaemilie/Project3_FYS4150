#ifndef __filename_hpp__
#define __filename_hpp__

#include <armadillo>

class Particle{
private:


public:
bool q;
bool m;
arma::vec r;
arma::vec v;

Particle(bool q, bool m, arma::vec r, arma::vec v){
    q = q;
    m = m;
    r = r;
    v = v;
};

};

int main()
{
    arma::vec r = arma::vec(1).fill(2.);
    arma::vec v = arma::vec(1).fill(2.);
    Particle p = Particle(1,1,r, v);
}

#endif