#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <vector>
#include <iostream>

//throw error if dx /= dy
//
//i = y
//j = x
//see https://isocpp.org/wiki/faq/templates
template <typename T> class Domain;
template <typename T> std::ostream& operator<< (std::ostream&, const std::vector<T>&);

template <class T> class Domain {
  public:
    Domain(); 
    ~Domain();
    Domain(int rows, int cols, int GC, double Lx); 

    int getNi(){return Ni;}; //computational rows
    int getNj(){return Nj;};
    std::vector<std::vector<T> > returnDomain() {return dom;};
    
    T operator() (unsigned int i, unsigned int j) const;

    //return row/column + ghost cells
    std::vector<T> getGRowi(int) const; // get the row in the ith direction
    std::vector<T> getGColj(int) const; //
    //return row/column
    std::vector<T> getRowi(int) const; //row 
    std::vector<T> getColj(int) const; //col
    int GCs(){return GC;};

    int starti, startj, endi, endj;
    int GNi, GNj;
    T getCellCartPos(int, int); //return the physical location of point indexed by cell centre i and j

    friend std::ostream& operator<< <>(std::ostream&, const std::vector<T>&);

  private:
    std::vector<std::vector<T> > dom;
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


