#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <vector>
#include <iostream>

//throw error if dx /= dy
//
//i = y
//j = x
class Domain {
  public:
    Domain(); 
    ~Domain();
    Domain(int rows, int cols, int GC, double Lx); 

    int getNi(){return Ni;}; //computational rows
    int getNj(){return Nj;};
    std::vector<std::vector<double> > returnDomain() {return dom;};
    
    double operator() (unsigned int i, unsigned int j) const;

    //return row/column + ghost cells
    std::vector<double> getGRowi(int) const; // get the row in the ith direction
    std::vector<double> getGColj(int) const; //
    //return row/column
    std::vector<double> getRowi(int) const; //row 
    std::vector<double> getColj(int) const; //col
    int GCs(){return GC;};

    int starti, startj, endi, endj;
    int GNi, GNj;
    double getPoint(int, int); //return the point indexed by i and j

    friend std::ostream& operator<<(std::ostream&, const std::vector<double>&);

  private:
    std::vector<std::vector<double> > dom;
    int GC;
    int Ni, Nj;
    //physical domain lengths and widths
    double Ly, Lx, dx;
    

    //solve once then take the transpose, dirn for the type of flux computed 
    //declare operators
  
    //get row, get col?
    //assign row, assign col.
  
};

#endif


