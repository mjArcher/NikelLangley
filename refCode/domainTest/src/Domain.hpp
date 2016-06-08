#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <vector>
#include <iostream>

//throw error if dx /= dy
class Domain {
  public:
    Domain(); 
    Domain(int rows, int cols, int GC, double Lx); 

    //number of cells excluding ghost cells 
    int getRows(){return Ny;};
    int getCols(){return Nx;};
    std::vector<std::vector<double> > returnDomain() {return dom;};

    int starti, endi;
    double Lx, Ly, dx;
    

    double operator() (unsigned int row, unsigned int col) const;

  private:
    std::vector<std::vector<double> > dom;
    int GC;
    int Nx, Ny;

    //solve once then take the transpose, dirn for the type of flux computed 
    //declare operators
  
    //get row, get col?
    //assign row, assign col.
  
};

#endif


