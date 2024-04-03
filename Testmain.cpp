
#include <iostream>
#include <armadillo>
#include "AO.h"
#include "utils.h"
#include "CNDO.h"


using namespace std;


int main() {

    AO inputfile_ao("H2.txt");

    // cout << "Overlap Matrix for H2: " << endl;
    vector<BasisFunction> basis_set = inputfile_ao.basis_set;

    map<int, BasisFunction> basis_map;
    for (int i = 0; i < basis_set.size(); i++) {
        if (basis_set[i].AO_type.find("s") != std::string::npos) {
            basis_map[basis_set[i].atom_index] = basis_set[i];
        }
    }

    arma::mat S = overlap_matrix(basis_set);
    // S.print();

    // create a CNDO instance
    CNDO CNDO_input;


    CNDO_input.updateDensityMatrix(inputfile_ao, "totalDensity");
    arma::mat S_deriv = overlapMatrix_derivative(basis_set);

    // Print the combined matrix

    cout << "Suv_RA" << endl;
    S_deriv.print();

    // cout << "Gamma_RA" << endl;
    // arma::field<arma::vec> gamma_deriv = CNDO_input.createGammaDerivativeMat(basis_map);
    // printField(gamma_deriv);





    return 0;
}
