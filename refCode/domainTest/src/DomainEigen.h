#ifndef DOMAINEIGEN_H
#define DOMAINEIGEN_H

#include <math.h>
#include <iomanip>
#include "ElasticPrimState.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>

typedef Eigen::Matrix<ElasticPrimState, Eigen::Dynamic, Eigen::Dynamic> Domain;

//function prototypes
void initialise(Domain&);
int afunction();

#endif 
