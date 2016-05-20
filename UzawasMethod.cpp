#include<vector>
#include<iostream>
#include<cmath>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>
#include <armadillo>

using namespace std;
using namespace arma;
//using namespace boost::numeric::ublas;


/*
Uzawa's method: Minimize the quadratic functional (Av,v) - (b,v) with the constraints: Cv <= d.
Here A must be symmetric positive definite and has dim = n.
Note A is square since it is symmetric.
There are m constraints and hence C is an mxn matrix.


Lagrange is an mx1 vector of lagrange multipliers
p is a parameter used in the convergence algorithm. It must be within certain bounds
to ensure convergence of the algorithm.
Some of the mat variables below are actually vectors. We leave them as nx1 mats for simplicity.

REFERENCE:
Introduction to Numerical Linear Algebra and Optimization. Philippe G. Ciarlet.
*/

mat UzawasMethod(mat A, mat b, mat C, mat d, double p, int iterations)
{
    
    int dim = A.n_rows;
    int dualDim = C.n_rows;
    mat Lagrange(dualDim,1);
    mat solution(dim,1);
    
    for(int i=0;i<iterations;i++)
    {
        for(int j=0;j<dualDim;j++)
        {
            mat temp = C*solution - d;
            Lagrange(j,0) = max(Lagrange(j,0) + p*(temp)(j,0),0.0);
        }
        solution = solve(A,b - C.t()*Lagrange);
    }
    return solution;
}

/*EXAMPLE 
Minimize 2x_1^2 + x_2^2 with respect to x_2 >= x_1 + 1

mat a(2,2);
    a(0,0) = 2; a(0,1) = 0;
    a(1,0) = 0; a(1,1) = 1;
    
mat b(2,1); 
    b(0,0) = 0; 
    b(1,0) = 0;

mat c(1,2); 
    c(0,0) = 1; c(0,1) = -1;

mat d(1,1); 
    d(0,0) = -1;
    
    mat s = UzawasMethod(a,b,c,d,.1,40);
    cout << s;
    
One can check directly using Lagrange Multiplier techniques that 
the solution given by the code for this example is indeed correct. 
x_1 = -1/3,  x_2 = 2/3
*/
