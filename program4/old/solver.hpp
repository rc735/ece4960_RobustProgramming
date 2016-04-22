//
//  solver.hpp
//  Program4
//
//  Created by Matthew Magaldi on 4/20/16.
//  Copyright © 2016 Matthew Magaldi. All rights reserved.
//

#ifndef solver_hpp
#define solver_hpp

#include <stdio.h>
#include <vector>

#include "methods.hpp"
#include "define.h"

typedef std::vector< double > state_type;
typedef std::vector< std::vector<double> > matrix;

matrix solver(double x, double t_final, int methodChoice);

#endif /* solver_hpp */
