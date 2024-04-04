#pragma once
#include <armadillo>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include "utils.h"

using namespace std;
using namespace arma;

struct basisFunction {
    string name;
    vec center;
    vec quantumNumbers;
    vec exponents;
    vec coefficients;
    vec normalization;
};

struct atom {
    string atomName;
    vec coords;
    int atomicNumber;
    vector<basisFunction> basisFunctionList;
    int numBasisSets_atom;
    int valence;
};

class Molecule {
public:
    int numAtoms;
    vector<atom> atomList;
    int numElectrons;
    int numBasisSets;
    int countElectrons();
    vector<basisFunction> basisFunctionListAll;
    int p;
    int q;
    int multiplicity;
    double totalEnergy;
    // create a matrix to store final density values 
    mat P_alpha_final;
    mat P_beta_final;
    double nuclear_repulsion_energy;
    double electronic_energy;

    // quantum numbers for the basis functions 
    vec s_lmn = { 0,0,0 };
    vec px_lmn = { 1,0,0 };
    vec py_lmn = { 0,1,0 };
    vec pz_lmn = { 0,0,1 };

    // exponents for the basis functions 
    vec H_exp = { 3.42525091, 0.62391373, 0.16885540 };
    vec C_exp = { 2.9412494, 0.6834831, 0.2222899 };
    vec N_exp = { 3.7804559, 0.8784966, 0.2857144 };
    vec O_exp = { 5.0331513, 1.1695961, 0.3803890 };
    vec F_exp = { 6.46480320, 1.50228120, 0.48858850 };

    // coefficients for the basis functions 
    vec H_coeff = { 0.15432897, 0.53532814, 0.44463454 };
    vec C_s_coeff = { -0.09996723, 0.39951283, 0.70011547 };
    vec C_p_coeff = { 0.15591627, 0.60768372, 0.39195739 };

    // Function: Build basis function for a given atom
    vector<basisFunction> buildBasisFunctions(atom& atom);

    // Constructor for the molecule class 
    Molecule(string filename);

    // print the geom info
    void print_geom_info();

    // Function: Calculate the normalization constants for a given basis function 
    void calculateNormalization(basisFunction& basis);

    // Function: Calculate the normalization constants for all of the basis functions
    void calcNormalizationAll(vector<basisFunction>& basisFunctionList);

    //Function: Calculate the contracted overlap integral for two primitive Gaussian shells
    double contractedOverlapIntegral(basisFunction basisFunction1, basisFunction basisFunction2);

    // Function: Calculate the overlap matrix of the entire molecule 
    mat overlapMatrix();

    // Function: Calculate the gamma matrix for two given basis functions 
    double calcGamma(basisFunction basisFunction1, basisFunction basisFunction2);

    // Function: Create the gamma matrix 
    mat calcGammaMatrix();

    // store the parameters needed for CNDO model 
    vec H_para = { 7.176, 0, 9 };
    vec C_para = { 14.051, 5.572, 21 };
    vec N_para = { 19.316, 7.275, 25 };
    vec O_para = { 25.390, 9.111, 31 };
    vec F_para = { 32.272, 11.080, 39 };

    // create a map to match the atomic number to the parameters
    map<string, vec> atom_para = { {"H", H_para}, {"C", C_para}, {"N", N_para}, {"O", O_para}, {"F", F_para} };

    // Function: Calulate the hamiltonian core matrix 
    mat calcHCoreMatrix();

    // Function: Update the density matrix 
    mat updateDensityMatrix(mat coeff, string spin);

    // Function : Calculate the total density for a certain atom using the P_alpha and P_beta matrices 
    vec totalDensity(mat P_alpha, mat P_beta);

    // Function: Calculate the Fock matrix values 
    mat calcFockMatrix(mat density_matrix, vec totalDensityVec);

    // Function: Calculate the nuclear repulsion energy 
    double calcNuclearRepulsionEnergy();

    // Function: Calculate the electronic energy
    double calcElectronicEnergy(mat density_matrix, mat Fock_matrix);

    // Function: classic SCF algorithm
    void runSCF(double tolerance, int maxIterations, int situation);

    // Function: Solve the X Matrix
    mat solveXMatrix(mat P_alpha, mat P_beta);

    // Function: Solve the Y Matrix 
    mat solveYMatrix(mat P_alpha, mat P_beta);

    vec contractedOverlapIntegral_derivative(basisFunction basisFunction1, basisFunction basisFunction2);

    // Function: Solve Overlap Integral RA 
    arma::field<arma::vec> overlapMatrix_derivative();

    double calc_derivative_vals(int dim, basisFunction basisFunction1, basisFunction basisFunction2);

    // Function: Calculate the 0 to 0 term 
    vec calc_zero_to_zero(basisFunction basisFunction1, basisFunction basisFunction2, double first_term, double second_term);

    // Function: Calculate the derivative of the elements in the gamma matrix 
    vec calcGammaDerivative(basisFunction basisFunction1, basisFunction basisFunction2);

    // Function: Calculate the derivative of the gamma matrix
    field<vec> calcGammaDerivativeMatrix();

    // Function: Calculate the gradient repulsion energy
    mat calcGradientRepulsionEnergy();

    // Function: Calculate the final gradient using the three values 
    mat calcFinalGradient();
};