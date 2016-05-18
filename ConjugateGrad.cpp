#include<vector>
#include<iostream>

using namespace std;

//Not complete will not work yet don't use

//Minimize the ellipic quadratic function J(v) = (Av,v) - (b,v) by conjugate gradient method with initial value u_0
//u_km1 stands for u_(k-1), the k-1(th) iteration
//w_k is the gradient of J(v) at u_k which has the alternate form written in the variable below
//d is the direction
//The while loop will terminate in finite time since this method is guaranteed to converge in n steps where n is the dimension

vector<double> conjugateGrad(vector<vector<double>> A, vector<double> b, vector<double> u_0)
{
    vector<double> d = times(A,u_0) - b //write matrix class to accomodate this
    double r = dot(d,d)/dot(times(A,d),d); //overload * operator for dot in mat class
    
    vector<double> u_k = u_0 - r*d;   // '-' or '*' not supported yet
    vector<double> u_km1 = u_0
    vector<double> w_k;
    vector<double> w_km1;
    
    while(times(A,u_k) - b != 0)
    {
        w_k =   A*u_k-b;
        w_km1 = A*u_km1-b;
        
        
        d = w_k + (w_k*w_k / w_km1*w_km1)*d;
        r = w_k*d / (A*d)*d;
        u_km1 = u_k;
        u_k = u_km1 - r*d;
    }
    return u_k;
}
