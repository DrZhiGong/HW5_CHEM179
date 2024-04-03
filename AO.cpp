
#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "AO.h"
#include "utils.h"

using namespace std;



//constructor
AO::AO(const char* filename) {

    //cout << "import AO object" << endl;
    ifstream infile(filename);
    vector<BasisFunction> basis_set;
    if (!infile.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        // Add additional error handling if needed
    }
    else {
        // Proceed with reading the file

        infile >> natoms;
        // make the next number the ionization 
        infile >> ionization;
        atom_types.resize(natoms);
        coord.resize(natoms, 3);  // this is an armadiilo matrix
        int atom_index = 0;
        for (int i = 0; i < natoms; i++) {
            atom_index = i;
            infile >> atom_types[i] >> coord(i, 0) >> coord(i, 1) >> coord(i, 2);
            // cout << "Atomic number: " << atom_types[i] << endl;
            if (atom_types[i] == "H") {
                num_hydrogen += 1;
                BasisFunction hydrogen_basis = createHydrogenBasis(atom_index, coord.row(i));
                // cout << "Hydrogen basis: " << hydrogen_basis.center << endl;
                basis_set.push_back(hydrogen_basis);
            }
            else if (atom_types[i] == "C") {
                num_carbon += 1;
                vector<BasisFunction> carbon_basis = createCarbonBasis(atom_index, coord.row(i));
                basis_set.insert(basis_set.end(), carbon_basis.begin(), carbon_basis.end());
            }
            else if (atom_types[i] == "N") {
                num_nitrogen += 1;
                vector<BasisFunction> nitrogen_basis = createNitrogenBasis(atom_index, coord.row(i));
                basis_set.insert(basis_set.end(), nitrogen_basis.begin(), nitrogen_basis.end());
            }
            else if (atom_types[i] == "O") {
                num_oxygen += 1;
                vector<BasisFunction> oxygen_basis = createOxygenBasis(atom_index, coord.row(i));
                basis_set.insert(basis_set.end(), oxygen_basis.begin(), oxygen_basis.end());
            }
            else if (atom_types[i] == "F") {
                num_fluorine += 1;
                vector<BasisFunction> fluorine_basis = createFluorineBasis(atom_index, coord.row(i));
                basis_set.insert(basis_set.end(), fluorine_basis.begin(), fluorine_basis.end());
            }
        }
        num_basis = count_basis(num_carbon, num_hydrogen, num_nitrogen, num_oxygen, num_fluorine);
        num_electrons = count_electrons(num_carbon, num_hydrogen, num_nitrogen, num_oxygen, num_fluorine, ionization);
    }

    updateBasisSet(basis_set);
    this->basis_set = basis_set;
    if (num_electrons % 2 == 0) {
        p = num_electrons / 2;
        q = num_electrons / 2;
    }
    else {
        p = (num_electrons + 1) / 2;
        q = (num_electrons - 1) / 2;
    }


    infile.close();
}

AO::~AO() {
    // destructor
}

int AO::count_basis(int num_carbon, int num_hydrogen, int num_nitrogen, int num_oxygen, int num_fluorine) {
    int num_basis = 0;
    num_basis = 4 * num_carbon + num_hydrogen + 4 * num_nitrogen + 4 * num_oxygen + 4 * num_fluorine;
    return num_basis;
}

int AO::count_electrons(int num_carbon, int num_hydrogen, int num_nitrogen, int num_oxygen, int num_fluorine, int ionization) {

    // returning a total number of valence electrons 
    int num_electrons = 4 * num_carbon + num_hydrogen + 5 * num_nitrogen + 6 * num_oxygen + 7 * num_fluorine - ionization;

    return num_electrons;
}

int AO::get_natoms() {
    return natoms;
}

arma::mat AO::get_coord() {
    return coord;
}

int AO::set_p(int new_p) {
    p = new_p;
    return p;
}


int AO::set_q(int new_q) {
    q = new_q;
    return q;
}


vector<string> AO::get_atom_types() {
    return atom_types;
}

// void AO::updateDensityMatrices(const arma::mat& new_density_alpha, const arma::mat& new_density_beta) {
//     density_alpha = new_density_alpha;
//     density_beta = new_density_beta;
// }

// void AO::updateCoefficientMatrices(const arma::mat& new_coeff_alpha, const arma::mat& new_coeff_beta) {
//     coefficient_alpha = new_coeff_alpha;
//     coefficient_beta = new_coeff_beta;
// }

// void AO::updateTotalVectors(const arma::vec& new_Ptotal) {
//     Ptotal = new_Ptotal;
// }

int AO::get_p() {
    return p;
}

int AO::get_q() {
    return q;
}


BasisFunction AO::createHydrogenBasis(int atom_index, arma::rowvec center) {
    return { "H_s", 1, {3.42525091, 0.62391373, 0.16885540}, {0.15432897, 0.53532814, 0.44463454}, center, {0, 0, 0}, {0.0, 0.0, 0.0}, atom_index };
}

vector<BasisFunction> AO::createCarbonBasis(int atom_index, arma::rowvec center) {
    vector<BasisFunction> carbonBasis;
    carbonBasis.push_back({ "C_s", 4, {2.94124940, 0.68348310, 0.22228990}, {-0.09996723, 0.39951283, 0.70011547}, center, {0, 0, 0}, {0.0, 0.0, 0.0}, atom_index });
    carbonBasis.push_back({ "C_px", 4, {2.94124940, 0.68348310, 0.22228990}, {0.15591627, 0.60768372, 0.39195739}, center, {1, 0, 0}, {0.0, 0.0, 0.0}, atom_index });
    carbonBasis.push_back({ "C_py", 4,  {2.94124940, 0.68348310, 0.22228990}, {0.15591627, 0.60768372, 0.39195739}, center, {0, 1, 0}, {0.0, 0.0, 0.0}, atom_index });
    carbonBasis.push_back({ "C_pz", 4, {2.94124940, 0.68348310, 0.22228990}, {0.15591627, 0.60768372, 0.39195739}, center, {0, 0, 1}, {0.0, 0.0, 0.0}, atom_index });
    return carbonBasis;
}

vector<BasisFunction> AO::createNitrogenBasis(int atom_index, arma::rowvec center) {
    vector<BasisFunction> nitrogenBasis;
    nitrogenBasis.push_back({ "N_s", 5, {3.78045588, 0.8784966, 0.2857144}, {-0.09996723, 0.39951283, 0.70011547}, center, {0, 0, 0}, {0.0, 0.0, 0.0}, atom_index });
    nitrogenBasis.push_back({ "N_px", 5, {3.78045588, 0.8784966, 0.2857144}, {0.15591627, 0.60768372, 0.39195739}, center, {1, 0, 0}, {0.0, 0.0, 0.0}, atom_index });
    nitrogenBasis.push_back({ "N_py", 5, {3.78045588, 0.8784966, 0.2857144}, {0.15591627, 0.60768372, 0.39195739}, center, {0, 1, 0}, {0.0, 0.0, 0.0}, atom_index });
    nitrogenBasis.push_back({ "N_pz", 5, {3.78045588, 0.8784966, 0.2857144}, {0.15591627, 0.60768372, 0.39195739}, center, {0, 0, 1}, {0.0, 0.0, 0.0}, atom_index });
    return nitrogenBasis;
}

vector<BasisFunction> AO::createOxygenBasis(int atom_index, arma::rowvec center) {
    vector<BasisFunction> oxygenBasis;
    oxygenBasis.push_back({ "O_s", 6, {5.03315130, 1.16959610, 0.38038900}, {-0.09996723, 0.39951283, 0.70011547}, center, {0, 0, 0}, {0.0, 0.0, 0.0}, atom_index });
    oxygenBasis.push_back({ "O_px", 6, {5.03315130, 1.16959610, 0.38038900}, {0.15591627, 0.60768372, 0.39195739}, center, {1, 0, 0}, {0.0, 0.0, 0.0}, atom_index });
    oxygenBasis.push_back({ "O_py", 6, {5.03315130, 1.16959610, 0.38038900}, {0.15591627, 0.60768372, 0.39195739}, center, {0, 1, 0}, {0.0, 0.0, 0.0}, atom_index });
    oxygenBasis.push_back({ "O_pz", 6, {5.03315130, 1.16959610, 0.38038900}, {0.15591627, 0.60768372, 0.39195739}, center, {0, 0, 1}, {0.0, 0.0, 0.0}, atom_index });
    return oxygenBasis;
}

vector<BasisFunction> AO::createFluorineBasis(int atom_index, arma::rowvec center) {
    vector<BasisFunction> fluorineBasis;
    fluorineBasis.push_back({ "F_s", 7, {6.46480320, 1.50228120, 0.48858850}, {-0.09996723, 0.39951283, 0.70011547}, center, {0, 0, 0}, {0.0, 0.0, 0.0}, atom_index });
    fluorineBasis.push_back({ "F_px", 7, {6.46480320, 1.50228120, 0.48858850}, {0.15591627, 0.60768372, 0.39195739}, center, {1, 0, 0}, {0.0, 0.0, 0.0}, atom_index });
    fluorineBasis.push_back({ "F_py", 7, {6.46480320, 1.50228120, 0.48858850}, {0.15591627, 0.60768372, 0.39195739}, center, {0, 1, 0}, {0.0, 0.0, 0.0}, atom_index });
    fluorineBasis.push_back({ "F_pz", 7, {6.46480320, 1.50228120, 0.48858850}, {0.15591627, 0.60768372, 0.39195739}, center, {0, 0, 1}, {0.0, 0.0, 0.0}, atom_index });
    return fluorineBasis;
}



void AO::updateNormalizationConstants(BasisFunction& basisFunction) {

    double omega1 = overlap_integral_3D(basisFunction.center, basisFunction.center, basisFunction.exponents(0), basisFunction.exponents(0), basisFunction.lmn, basisFunction.lmn);
    double omega2 = overlap_integral_3D(basisFunction.center, basisFunction.center, basisFunction.exponents(1), basisFunction.exponents(1), basisFunction.lmn, basisFunction.lmn);
    double omega3 = overlap_integral_3D(basisFunction.center, basisFunction.center, basisFunction.exponents(2), basisFunction.exponents(2), basisFunction.lmn, basisFunction.lmn);

    // cout << "omega1: " << omega1 << endl;
    // cout << "omega2: " << omega2 << endl;
    // cout << "omega3: " << omega3 << endl;
    basisFunction.normalization_constants(0) = pow(omega1, -0.5);
    basisFunction.normalization_constants(1) = pow(omega2, -0.5);
    basisFunction.normalization_constants(2) = pow(omega3, -0.5);

    // cout << "normalization_constants: " << basisFunction.normalization_constants << endl;

}




void AO::updateBasisSet(vector<BasisFunction>& basis_set) {
    for (int i = 0; i < basis_set.size(); i++) {
        updateNormalizationConstants(basis_set[i]);
    }
}



double overlap_integral_1D(double center1, double center2, double alpha1, double alpha2, double l1, double l2) {

    double prefactor_1D = prefactor(center1, alpha1, center2, alpha2);
    double center_product = center_of_product(center1, alpha1, center2, alpha2);
    double double_sum = 0.0;

    for (int i = 0; i <= l1; i++) {
        for (int j = 0; j <= l2; j++) {
            if ((i + j) % 2 == 0) {
                double_sum += binomial(l1, i) * binomial(l2, j) * ((double_factorial(i + j - 1) * pow(center_product - center1, l1 - i) * pow(center_product - center2, l2 - j)) / pow(2 * (alpha1 + alpha2), double(i + j) / 2));
            }
        }
    }

    return double_sum * prefactor_1D;
}


double overlap_integral_3D(arma::rowvec center1, arma::rowvec center2, double alpha1, double alpha2, arma::vec lmn1, arma::vec lmn2) {

    double overlap_orbital = overlap_integral_1D(center1(0), center2(0), alpha1, alpha2, lmn1(0), lmn2(0)) \
        * overlap_integral_1D(center1(1), center2(1), alpha1, alpha2, lmn1(1), lmn2(1)) \
        * overlap_integral_1D(center1(2), center2(2), alpha1, alpha2, lmn1(2), lmn2(2));

    // cout << "used in normalization: " << overlap_orbital << endl;
    return overlap_orbital;

}
double unnormOverlapIntegral(BasisFunction& basisFunction1, BasisFunction& basisFunction2, int i, int j) {

    double overlap_orbital = overlap_integral_1D(basisFunction1.center(0), basisFunction2.center(0), basisFunction1.exponents(i), basisFunction2.exponents(j), basisFunction1.lmn(0), basisFunction2.lmn(0)) \
        * overlap_integral_1D(basisFunction1.center(1), basisFunction2.center(1), basisFunction1.exponents(i), basisFunction2.exponents(j), basisFunction1.lmn(1), basisFunction2.lmn(1)) \
        * overlap_integral_1D(basisFunction1.center(2), basisFunction2.center(2), basisFunction1.exponents(i), basisFunction2.exponents(j), basisFunction1.lmn(2), basisFunction2.lmn(2));
    // cout << "Unnormalized overlap_orbital: " << overlap_orbital << endl;
    return overlap_orbital;
}

double overlapIntegral(BasisFunction& basisFunction1, BasisFunction& basisFunction2) {

    double overlap_orbital = 0.0;
    // double unnormalized_overlap = unnormOverlapIntegral(basisFunction1, basisFunction2);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {

            overlap_orbital += basisFunction1.normalization_constants(i) * basisFunction2.normalization_constants(j) *
                basisFunction1.contracted_coefficients(i) * basisFunction2.contracted_coefficients(j) * unnormOverlapIntegral(basisFunction1, basisFunction2, i, j);
        }
    }


    // cout << "Normalized overlap_orbital: " << overlap_orbital << endl;

    return overlap_orbital;
}

arma::mat overlap_matrix(vector<BasisFunction>& basis_set) {
    arma::mat S = arma::zeros(basis_set.size(), basis_set.size());
    for (int i = 0; i < basis_set.size(); i++) {
        for (int j = 0; j < basis_set.size(); j++) {

            double temp;
            if (abs(overlapIntegral(basis_set[i], basis_set[j])) < 1e-15) {
                temp = 0.0;
            }
            else {
                temp = overlapIntegral(basis_set[i], basis_set[j]);
            }

            S(i, j) = temp;
            // cout << "inside overlap matrix call overlap integral. DONE!" << endl;
        }
    }
    return S;
}


double derivative_OV_1D(double center1, double center2, double alpha1, double alpha2, int l1, int l2) {


    double derivative = 0.0;

    if (l1 == 0) {
        derivative = 2 * alpha1 * overlap_integral_1D(center1, center2, alpha1, alpha2, l1 + 1, l2);
    }

    if (l1 == 1) {
        derivative = -l1 * overlap_integral_1D(center1, center2, alpha1, alpha2, l1 - 1, l2) +
            2 * alpha1 * overlap_integral_1D(center1, center2, alpha1, alpha2, l1 + 1, l2);
    }

    return derivative;

}

arma::vec derivative_OV_3D(arma::rowvec centers_a, arma::rowvec centers_b, double alpha, double alpha1, arma::vec lmnA, arma::vec lmnB) {

    double x = overlap_integral_1D(centers_a(0), centers_b(0), alpha, alpha1, lmnA(0), lmnB(0));
    double y = overlap_integral_1D(centers_a(1), centers_b(1), alpha, alpha1, lmnA(1), lmnB(1));
    double z = overlap_integral_1D(centers_a(2), centers_b(2), alpha, alpha1, lmnA(2), lmnB(2));

    // Calculate the overlap integral gradient for each dimension
    double dx = derivative_OV_1D(centers_a(0), centers_b(0), alpha, alpha1, lmnA(0), lmnB(0));
    double dy = derivative_OV_1D(centers_a(0), centers_b(0), alpha, alpha1, lmnA(1), lmnB(1));
    double dz = derivative_OV_1D(centers_a(0), centers_b(0), alpha, alpha1, lmnA(2), lmnB(2));

    arma::vec ov_gradient = { dx * y * z, x * dy * z, x * y * dz };
    return ov_gradient;
}

arma::vec contractedOV_derivative(BasisFunction basisFunction1, BasisFunction basisFunction2) {
    arma::vec derivatives = arma::zeros(3);

    for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
            derivatives += basisFunction1.normalization_constants(k) * basisFunction2.normalization_constants(l) *
                basisFunction1.contracted_coefficients(k) * basisFunction2.contracted_coefficients(l) *
                derivative_OV_3D(basisFunction1.center, basisFunction2.center, basisFunction1.exponents(k), basisFunction2.exponents(l), basisFunction1.lmn, basisFunction2.lmn);
        }
    }

    return derivatives;
}

arma::mat overlapMatrix_derivative(vector<BasisFunction>& basis_set) {

    // std::cout << "inside overlapMatrix_derivative" << std::endl;
    arma::mat overlapDerivativeMat = arma::zeros(3, basis_set.size() * basis_set.size());

    for (int mu = 0; mu < basis_set.size(); mu++) {
        for (int nu = 0; nu < basis_set.size(); nu++) {

            if (basis_set[mu].atom_index != basis_set[nu].atom_index) {
                // std::cout << "inside overlapMatrix_derivative if statement" << std::endl;
                overlapDerivativeMat.col(mu * basis_set.size() + nu) = contractedOV_derivative(basis_set[mu], basis_set[nu]);
                // std::cout << "result: " << result << std::endl;
            }
            else {
                overlapDerivativeMat.col(mu * basis_set.size() + nu) = arma::zeros(3);
            }


        }
    }
    return overlapDerivativeMat;
}















