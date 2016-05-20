#include<vector>
#include<iostream>
#include<cmath>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>
#include <armadillo>

using namespace std;
using namespace arma;
//using namespace boost::numeric::ublas;


//Uzawa's method. Let A be symmetric positive definite with dim = n.
//Note A is square since it is symmetric.
//Here C is an mxn matrix of constraints
//Lagrange is an mx1 vector of lagrange multipliers
//p is a parameter used in the convergence algorithm. It must be within certain bounds
//to ensure convergence of the algorithm.
//Some of the mat variables below are actually vectors. 
//We leave them as nx1 mats for simplicity.

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

mat a(2,2);
    a(0,0) = 2; a(0,1) = 0;
    a(1,0) = 0; a(1,1) = 1;
    mat b(2,1); b(0,0) = 0; b(1,0) = 0;
    mat c(1,2); c(0,0) = 1; c(0,1) = -1;
    mat d(1,1); d(0,0) = -1;
    
    mat s = UzawasMethod(a,b,c,d,.1,40);
    cout << s;
    
One can check directly using simple Lagrange Multiplier techniques that 
the solution given by the code for this example is indeed correct. 

*/
