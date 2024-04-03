#pragma once

#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>


using namespace std;

struct BasisFunction {
    string AO_type;
    int valence_e;
    arma::vec exponents;
    arma::vec contracted_coefficients;
    arma::rowvec center;
    arma::vec lmn;
    arma::vec normalization_constants;
    int atom_index;
};

class AO {

private:
    int num_hydrogen;
    int num_carbon;
    int num_nitrogen;
    int num_oxygen;
    int num_fluorine;
    int num_basis;
    int natoms;
    int ionization;
    int p;
    int q;
    arma::mat coord;
    vector<string> atom_types;


    int count_basis(int num_carbon, int num_hydrogen, int num_nitrogen, int num_oxygen, int num_fluorine);  //count the number of basis functions
    int count_electrons(int num_carbon, int num_hydrogen, int num_nitrogen, int num_oxygen, int num_fluorine, int ionization); //count the number of electrons
    BasisFunction createHydrogenBasis(int atom_index, arma::rowvec center); //create a hydrogen basis function
    vector<BasisFunction> createCarbonBasis(int atom_index, arma::rowvec center); //create a carbon basis function
    vector<BasisFunction> createNitrogenBasis(int atom_index, arma::rowvec center); //create a nitrogen basis function
    vector<BasisFunction> createOxygenBasis(int atom_index, arma::rowvec center); //create a oxygen basis function
    vector<BasisFunction> createFluorineBasis(int atom_index, arma::rowvec center); //create a fluorine basis function
    void updateNormalizationConstants(BasisFunction& basisFunction);
    void updateBasisSet(vector<BasisFunction>& basis_set);



public:
    //constructor
    AO(const char* filename);
    ~AO();
    int num_electrons;
    // arma::mat density_alpha;
    // arma::mat density_beta;
    // arma::mat coefficient_alpha;
    // arma::mat coefficient_beta;
    // arma::vec Ptotal;
    int get_natoms();
    vector<string> get_atom_types();

    int get_p();
    int get_q();
    arma::mat get_coord();
    int set_p(int new_p);
    int set_q(int new_q);

    vector<BasisFunction> basis_set;


};

double overlap_integral_1D(double center1, double center2, double alpha1, double alpha2, double l1, double l2);
double overlap_integral_3D(arma::rowvec center1, arma::rowvec center2, double alpha1, double alpha2, arma::vec lmn1, arma::vec lmn2);
double unnormOverlapIntegral(BasisFunction& basisFunction1, BasisFunction& basisFunction2, int i, int j);
double overlapIntegral(BasisFunction& basisFunction1, BasisFunction& basisFunction2);
arma::mat overlap_matrix(vector<BasisFunction>& basis_set);
double derivative_OV_1D(double center1, double center2, double alpha1, double alpha2, int l1, int l2);
arma::vec derivative_OV_3D(arma::rowvec centers_a, arma::rowvec centers_b, double alpha, double alpha1, arma::vec lmnA, arma::vec lmnB);
arma::vec contractedOV_derivative(BasisFunction basisFunction1, BasisFunction basisFunction2);
arma::mat overlapMatrix_derivative(vector<BasisFunction>& basis_set);

