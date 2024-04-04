#pragma once 
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

// Function: Evaluate the double factorial of a number
int double_factorial(int n);

// Function: Evaluate the exponential prefactor of the product of the Gaussian
double product_exponential_prefactor(double xa, double xb, double alpha1, double alpha2);

// Function: Evaluate the factorial of a given number 
int factorial(int n);

// Function: Evaluate the binomial prefactor of the product of the Gaussian
int product_binomial_prefactor(int i, int x);

// Function: Evaluate the overlap integral at a specific dimension 
double overlap_at_1D(double alpha1, double alpha2, double xa, double xb, int l1, int l2);

//Function: Evaluate the overlap integral for two primitive Gaussian shells
double overlapIntegral(vec center1, vec center2, double alpha1, double alpha2, vec lmn1, vec lmn2);

//Function: Evaluate the overlap integral derivative for two primitive Gaussian shells
double overlapIntegral_derivative(vec center1, vec center2, double alpha1, double alpha2, vec lmn1, vec lmn2, int coordinate);

// Function: Evaluate the overlap integral at a specific dimension 
double overlap_at_1D(double alpha1, double alpha2, double xa, double xb, int l1, int l2, int coordinate);

// Function: Calculate the derivative factor for the specified coordinate
double gradient_factor(int i, int j, int coordinate, double xa, double xb, double xp, int l1, int l2);

// Function: Helper function calculate the 2e integral at a specific dimension
double I2e_pG(vec& Ra, vec Rb, double sigmaA, double sigmaB);

// Function: Swap coordinates of a vector 
void swapCoordinates(vec& v, int dimension);

//Function: Calculate the 2e integral derivative at a specific dimension
vec I2e_pG_derivative(vec& Ra, vec Rb, double sigmaA, double sigmaB);

void print_field(arma::field<arma::vec> field);