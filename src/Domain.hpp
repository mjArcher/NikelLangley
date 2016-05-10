#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <vector>
#include <iostream>

class Domain {
  public:
    Domain(); 
    Domain(int rows, int cols); 
    int getRows(){return rows;};
    int getCols(){return cols;};
    // this is bad but i'm doing it:
    std::vector<std::vector<double> > returnDomain(){return dom;};

    double operator() (unsigned int row, unsigned int col) const;


  private:
    std::vector<std::vector<double> > dom;
    int rows;
    int cols;
  //declare operators
  
    //get row, get col?
    //assign row, assign col.

  
};

#endif
