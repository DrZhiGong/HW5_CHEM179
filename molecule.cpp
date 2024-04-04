#include "molecule.h"
#include <fstream> 

const double M_PI = 3.1415926;
const double Hartree = 27.211396641308;

// Constructor for the molecule class 
Molecule::Molecule(string filename) {
    ifstream infile(filename);
    if (!infile) {
        cout << "File does not exist" << endl;
        exit(1);
    }

    // map for the elements 
    map<string, int> element = { {"H", 1}, {"C",6}, {"N",7}, {"O",8}, {"F",9} };
    map<string, int> element_valence = { {"H", 1}, {"C",4}, {"N",5}, {"O",6}, {"F",7} };


    infile >> numAtoms >> multiplicity;
    atomList = vector<atom>(numAtoms);
    for (int i = 0; i < numAtoms; i++) {
        // initialize the atom vector coord 
        atomList[i].coords = zeros<vec>(3);
        infile >> atomList[i].atomName >> atomList[i].coords(0) >> atomList[i].coords(1) >> atomList[i].coords(2);
        atomList[i].atomicNumber = element[atomList[i].atomName];
        atomList[i].valence = element_valence[atomList[i].atomName];
    }


    infile.close();

    // count the number of electrons 
    numElectrons = countElectrons();

    // build basis functions for all the atoms 
    numBasisSets = 0;
    for (int i = 0; i < numAtoms; i++) {
        atomList[i].basisFunctionList = buildBasisFunctions(atomList[i]);
        // calculate the normalization constant
        calcNormalizationAll(atomList[i].basisFunctionList);
        numBasisSets += atomList[i].basisFunctionList.size();
        for (int j = 0; j < atomList[i].basisFunctionList.size(); j++) {
            basisFunctionListAll.push_back(atomList[i].basisFunctionList[j]);
        }
    }

}

// Function: Count the number of electrons in the molecule 
int Molecule::countElectrons() {
    int count = 0;
    for (int i = 0; i < numAtoms; i++) {
        if (atomList[i].atomicNumber == 1) {
            count += 1;
        }
        else if (atomList[i].atomicNumber == 6) {
            count += 4;
        }
        else if (atomList[i].atomicNumber == 7) {
            count += 5;
        }
        else if (atomList[i].atomicNumber == 8) {
            count += 6;
        }
        else if (atomList[i].atomicNumber == 9) {
            count += 7;
        }
    }

    // multiplicity 
    if (multiplicity == 0) {
        p = count / 2;
        q = count - p;
    }
    if (multiplicity == 1) {
        count -= 1;
        p = count / 2;
        q = count - p;
    }
    if (multiplicity == 3) {
        p = 4;
        q = 1;
    }



    return count / 2;
}

// Build basis functions 
vector<basisFunction> Molecule::buildBasisFunctions(atom& atom) {
    vector<basisFunction> basisFunctionList;

    int atomNum = atom.atomicNumber;
    if (atomNum == 1) {
        struct basisFunction s = { "1s", atom.coords, s_lmn, H_exp, H_coeff, zeros<vec>(3) };
        basisFunctionList.push_back(s);
    }
    else if (atomNum == 6) {
        struct basisFunction twos = { "2s", atom.coords, s_lmn, C_exp, C_s_coeff, zeros<vec>(3) };
        struct basisFunction px = { "2px", atom.coords, px_lmn, C_exp, C_p_coeff, zeros<vec>(3) };
        struct basisFunction py = { "2py", atom.coords, py_lmn, C_exp, C_p_coeff, zeros<vec>(3) };
        struct basisFunction pz = { "2pz", atom.coords, pz_lmn, C_exp, C_p_coeff, zeros<vec>(3) };
        basisFunctionList.push_back(twos);
        basisFunctionList.push_back(px);
        basisFunctionList.push_back(py);
        basisFunctionList.push_back(pz);
    }
    else if (atomNum == 7) {
        struct basisFunction twos = { "2s", atom.coords, s_lmn, N_exp, C_s_coeff, zeros<vec>(3) };
        struct basisFunction px = { "2px", atom.coords, px_lmn, N_exp, C_p_coeff, zeros<vec>(3) };
        struct basisFunction py = { "2py", atom.coords, py_lmn, N_exp, C_p_coeff, zeros<vec>(3) };
        struct basisFunction pz = { "2pz", atom.coords, pz_lmn, N_exp, C_p_coeff, zeros<vec>(3) };
        basisFunctionList.push_back(twos);
        basisFunctionList.push_back(px);
        basisFunctionList.push_back(py);
        basisFunctionList.push_back(pz);
    }
    else if (atomNum == 8) {
        struct basisFunction twos = { "2s", atom.coords, s_lmn, O_exp, C_s_coeff, zeros<vec>(3) };
        struct basisFunction px = { "2px", atom.coords, px_lmn, O_exp, C_p_coeff, zeros<vec>(3) };
        struct basisFunction py = { "2py", atom.coords, py_lmn, O_exp, C_p_coeff, zeros<vec>(3) };
        struct basisFunction pz = { "2pz", atom.coords, pz_lmn, O_exp, C_p_coeff, zeros<vec>(3) };
        basisFunctionList.push_back(twos);
        basisFunctionList.push_back(px);
        basisFunctionList.push_back(py);
        basisFunctionList.push_back(pz);
    }
    else if (atomNum == 9) {
        struct basisFunction twos = { "2s", atom.coords, s_lmn, F_exp, C_s_coeff, zeros<vec>(3) };
        struct basisFunction px = { "2px", atom.coords, px_lmn, F_exp, C_p_coeff, zeros<vec>(3) };
        struct basisFunction py = { "2py", atom.coords, py_lmn, F_exp, C_p_coeff, zeros<vec>(3) };
        struct basisFunction pz = { "2pz", atom.coords, pz_lmn, F_exp, C_p_coeff, zeros<vec>(3) };
        basisFunctionList.push_back(twos);
        basisFunctionList.push_back(px);
        basisFunctionList.push_back(py);
        basisFunctionList.push_back(pz);
    }

    return basisFunctionList;
}



// Function: Calculate the Normalization Constants for a given basis function 
void Molecule::calculateNormalization(basisFunction& basisFunction) {
    // Calculate the normalization constant for the basis function
    double norm1 = overlapIntegral(basisFunction.center, basisFunction.center, basisFunction.exponents(0), basisFunction.exponents(0), basisFunction.quantumNumbers, basisFunction.quantumNumbers);
    double norm2 = overlapIntegral(basisFunction.center, basisFunction.center, basisFunction.exponents(1), basisFunction.exponents(1), basisFunction.quantumNumbers, basisFunction.quantumNumbers);
    double norm3 = overlapIntegral(basisFunction.center, basisFunction.center, basisFunction.exponents(2), basisFunction.exponents(2), basisFunction.quantumNumbers, basisFunction.quantumNumbers);

    norm1 = pow(norm1, -0.5);
    norm2 = pow(norm2, -0.5);
    norm3 = pow(norm3, -0.5);

    basisFunction.normalization = { norm1, norm2, norm3 };
}

// Function: Calculate the Normalization Constants for all of the basis functions
void Molecule::calcNormalizationAll(vector<basisFunction>& basisFunctionList) {
    for (int i = 0; i < basisFunctionList.size(); i++) {
        calculateNormalization(basisFunctionList[i]);
    }
}

//Function: Calculate the contracted overlap integral for two primitive Gaussian shells
double Molecule::contractedOverlapIntegral(basisFunction basisFunction1, basisFunction basisFunction2) {
    double overlap = 0;

    // loop over all the exponent combinations 
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            overlap += basisFunction1.coefficients(i) * basisFunction2.coefficients(j) * basisFunction1.normalization(i) * basisFunction2.normalization(j) *
                overlapIntegral(basisFunction1.center, basisFunction2.center, basisFunction1.exponents(i), basisFunction2.exponents(j), basisFunction1.quantumNumbers, basisFunction2.quantumNumbers);
        }
    }

    return overlap;
}

// Function: Calculate the overlap integral for two primitive Gaussian shells
mat Molecule::overlapMatrix() {
    mat overlap = arma::zeros<mat>(numBasisSets, numBasisSets);
    for (int i = 0; i < numBasisSets; i++) {
        for (int j = 0; j < numBasisSets; j++) {
            overlap(i, j) = contractedOverlapIntegral(basisFunctionListAll[i], basisFunctionListAll[j]);
        }
    }
    return overlap;
}

// Function: Calculate the gamma value for two given basis functions
double Molecule::calcGamma(basisFunction basisFunction1, basisFunction basisFunction2) {
    // get the lmn values 
    vec lmn1 = basisFunction1.quantumNumbers;
    vec lmn2 = basisFunction2.quantumNumbers;

    // make sure only the s orbitals are being used
    if (lmn1(0) != 0 || lmn1(1) != 0 || lmn1(2) != 0) {
        cout << "Error: Only s orbitals are allowed" << endl;
    }

    // extract other information from the basis sets
    vec da = basisFunction1.coefficients % basisFunction1.normalization;
    vec db = basisFunction2.coefficients % basisFunction2.normalization;
    vec alphaa = basisFunction1.exponents;
    vec alphab = basisFunction2.exponents;
    vec Ra = basisFunction1.center;
    vec Rb = basisFunction2.center;
    int len = basisFunction1.exponents.size();

    double sum = 0;
    for (int k1 = 0; k1 < len; k1++) {
        for (int k2 = 0; k2 < len; k2++) {
            double sigmaA = 1.0 / (alphaa(k1) + alphaa(k2));

            for (int j1 = 0; j1 < len; j1++) {
                for (int j2 = 0; j2 < len; j2++) {
                    double sigmaB = 1.0 / (alphab(j1) + alphab(j2));

                    double I2e = I2e_pG(Ra, Rb, sigmaA, sigmaB);

                    sum += da(k1) * da(k2) * db(j1) * db(j2) * I2e;
                }
            }
        }
    }

    return sum * Hartree;
}

// Function: create the gamma matrix 
mat Molecule::calcGammaMatrix() {
    mat gamma = arma::zeros<mat>(numAtoms, numAtoms);
    for (int i = 0; i < numAtoms; i++) {
        for (int j = 0; j < numAtoms; j++) {
            gamma(i, j) = calcGamma(atomList[i].basisFunctionList[0], atomList[j].basisFunctionList[0]);
        }
    }
    return gamma;
}

// Function: Calculate the hamiltonian core matrix 
mat Molecule::calcHCoreMatrix() {
    mat hcore = arma::zeros<mat>(numBasisSets, numBasisSets);
    int mu = 0;
    double para_val = 0;

    for (int A = 0; A < numAtoms; A++) {
        string name_A = atomList[A].atomName;
        int z_val_A = atomList[A].valence;
        int numBasis_A = atomList[A].basisFunctionList.size();
        double gammaAA = calcGamma(atomList[A].basisFunctionList[0], atomList[A].basisFunctionList[0]);
        // cout << "gammaAA: " << gammaAA << endl;

        // loop over all the basis sets in atomA
        for (int i = 0; i < numBasis_A; i++) {
            string name_ao = atomList[A].basisFunctionList[i].name;
            // find the related parameter 
            vec para = atom_para[name_A];
            if (name_ao == "1s" || name_ao == "2s") {
                para_val = para(0);
            }
            else if (name_ao == "2px" || name_ao == "2py" || name_ao == "2pz") {
                para_val = para(1);
            }
            hcore(mu, mu) = -para_val - (z_val_A - 0.5) * gammaAA;

            int nu = 0;

            // loop over all the atoms in the molecule 
            for (int B = 0; B < numAtoms; B++) {
                string name_B = atomList[B].atomName;
                int z_val_B = atomList[B].valence;
                int numBasis_B = atomList[B].basisFunctionList.size();
                double gammaAB = calcGamma(atomList[A].basisFunctionList[0], atomList[B].basisFunctionList[0]);

                // add the ZB and gammaAB terms
                if (A != B) {
                    hcore(mu, mu) -= (z_val_B * gammaAB);
                }

                // loop over all the basis sets in atomB 
                for (int j = 0; j < numBasis_B; j++) {
                    if (mu != nu) {
                        int beta_A = atom_para[name_A](2);
                        int beta_B = atom_para[name_B](2);
                        hcore(mu, nu) = -0.5 * (beta_A + beta_B) * contractedOverlapIntegral(atomList[A].basisFunctionList[i], atomList[B].basisFunctionList[j]);
                    }
                    nu += 1;
                }
            }
            mu += 1;
        }
    }
    return hcore;
}

//Function: Update the density matrix 
mat Molecule::updateDensityMatrix(mat coeff, string spin) {
    mat matrix = zeros<mat>(numBasisSets, numBasisSets);
    // p(mu,vu) = sum over p electrons of C(p,mu)*C(p,vu)
    if (spin == "alpha") {
        for (int mu = 0; mu < numBasisSets; mu++) {
            for (int vu = 0; vu < numBasisSets; vu++) {
                for (int i = 0; i < p; i++) {
                    matrix(mu, vu) += coeff(mu, i) * coeff(vu, i);
                }
            }
        }
        // p(mu,vu) = sum over q electrons of C(q,mu)*C(q,vu)
    }
    else if (spin == "beta") {
        for (int mu = 0; mu < numBasisSets; mu++) {
            for (int vu = 0; vu < numBasisSets; vu++) {
                for (int i = 0; i < q; i++) {
                    matrix(mu, vu) += coeff(mu, i) * coeff(vu, i);
                }
            }
        }
    }
    return matrix;
}

// Function: Calculate the total density for a certain atom using the P_alpha and P_beta matrices 
vec Molecule::totalDensity(mat P_alpha, mat P_beta) {
    vec totalDensity = zeros<vec>(numAtoms);

    // add alpha and beta density together
    mat newMatrix = P_alpha + P_beta;
    int index = 0;

    for (int i = 0; i < numAtoms; i++) {
        int numBasis_A = atomList[i].basisFunctionList.size();
        for (int j = 0; j < numBasis_A; j++) {
            totalDensity(i) += newMatrix(index, index);
            index += 1;
        }
    }

    return totalDensity;
}

// Function: Calculate the fock matric for alpha electrons
mat Molecule::calcFockMatrix(mat density_matrix, vec totalDensityVec) {
    mat fock = zeros<mat>(numBasisSets, numBasisSets);
    int mu = 0;
    double para_val = 0;

    for (int A = 0; A < numAtoms; A++) {
        string name_A = atomList[A].atomName;
        int z_val_A = atomList[A].valence;
        int numBasis_A = atomList[A].basisFunctionList.size();
        double gammaAA = calcGamma(atomList[A].basisFunctionList[0], atomList[A].basisFunctionList[0]);
        double density = totalDensityVec(A);


        // loop over all the basis sets in atomA
        for (int i = 0; i < numBasis_A; i++) {
            string name_ao = atomList[A].basisFunctionList[i].name;
            // find the related parameter 
            vec para = atom_para[name_A];
            if (name_ao == "1s" || name_ao == "2s") {
                para_val = para(0);
            }
            else if (name_ao == "2px" || name_ao == "2py" || name_ao == "2pz") {
                para_val = para(1);
            }

            fock(mu, mu) = -para_val + ((density - z_val_A) - (density_matrix(mu, mu) - 0.5)) * gammaAA;

            int nu = 0;

            // loop over all the atoms in the molecule 
            for (int B = 0; B < numAtoms; B++) {
                string name_B = atomList[B].atomName;
                int z_val_B = atomList[B].valence;
                int numBasis_B = atomList[B].basisFunctionList.size();
                double gammaAB = calcGamma(atomList[A].basisFunctionList[0], atomList[B].basisFunctionList[0]);
                double density = totalDensityVec(B);

                // add the ZB and gammaAB terms
                if (A != B) {
                    fock(mu, mu) += (density - z_val_B) * gammaAB;
                }

                // loop over all the basis sets in atomB 
                for (int j = 0; j < numBasis_B; j++) {
                    if (mu != nu) {
                        int beta_A = atom_para[name_A](2);
                        int beta_B = atom_para[name_B](2);
                        fock(mu, nu) = -0.5 * (beta_A + beta_B) * contractedOverlapIntegral(atomList[A].basisFunctionList[i], atomList[B].basisFunctionList[j])
                            - (density_matrix(mu, nu) * gammaAB);
                    }
                    nu += 1;
                }
            }
            mu += 1;
        }
    }
    return fock;
}

// Function: Calculate nuclear repulsion energy for the molecule
double Molecule::calcNuclearRepulsionEnergy() {
    double energy = 0;
    for (int i = 0; i < numAtoms; i++) {
        for (int j = i + 1; j < numAtoms; j++) {
            double distance = norm(atomList[i].coords - atomList[j].coords);
            energy += (atomList[i].valence * atomList[j].valence) / distance;
        }
    }
    return energy;
}

// Function: Calculate electronic repulsion 
double Molecule::calcElectronicEnergy(mat density_matrix, mat Fock_matrix) {
    double energy = 0;
    mat hcore = calcHCoreMatrix();
    for (int mu = 0; mu < numBasisSets; mu++) {
        for (int nu = 0; nu < numBasisSets; nu++) {
            energy += density_matrix(mu, nu) * (hcore(mu, nu) + Fock_matrix(mu, nu));
        }
    }
    return energy * 0.5;
}

// Function: Calculate the classic SCF algorithm 
void Molecule::runSCF(double tolerance, int maxIterations,int situation) {
    // step 1: Guess the density matrices 
    mat P_alpha = zeros<mat>(numBasisSets, numBasisSets);
    mat P_beta = zeros<mat>(numBasisSets, numBasisSets);
    vec totalDensityvec = zeros<vec>(numAtoms);
    atomList[1].coords(0) = atomList[1].coords(0) + situation * 1e-6;

    bool converged = false;
    int iteration = 0;


    // step 6: Check if the density matrices have converged
    while (!converged && iteration < maxIterations) {

        //cout << "Iteration: " << iteration << endl;

        // step 2: Calculate the Fock matrix 
        mat F_alpha = calcFockMatrix(P_alpha, totalDensityvec);
        mat F_beta = calcFockMatrix(P_beta, totalDensityvec);

        // step 3: diagonalize the Fock matrix and obtain the eigenvalues and eigenvectors
        vec eigenvalues_alpha;
        mat eigenvectors_alpha;
        eig_sym(eigenvalues_alpha, eigenvectors_alpha, F_alpha);

        vec eigenvalues_beta;
        mat eigenvectors_beta;
        eig_sym(eigenvalues_beta, eigenvectors_beta, F_beta);

        //cout << "Solving the eigenvalue problem : " << iteration << endl;


        // save the old density matrices
        mat P_alpha_old = P_alpha;
        mat P_beta_old = P_beta;



        // step 4: Calculate the new density matrices
        P_alpha = updateDensityMatrix(eigenvectors_alpha, "alpha");


        P_beta = updateDensityMatrix(eigenvectors_beta, "beta");


        totalDensityvec = totalDensity(P_alpha, P_beta);


        // step 5: check for convergence 
        if (abs(P_alpha - P_alpha_old).max() < tolerance && abs(P_beta - P_beta_old).max() < tolerance) {
            converged = true;
            if (converged) {

                cout << "Nuclear Repulsion Energy is " ;
                double nre = calcNuclearRepulsionEnergy() * Hartree;
                cout << nre << " eV" << endl;
                double ee_alpha = calcElectronicEnergy(P_alpha, F_alpha);
                double ee_beta = calcElectronicEnergy(P_beta, F_beta);
                cout << "Electron Energy is";
               cout << (ee_alpha + ee_beta) << " eV" << endl;
                cout << "(Not required for HW5)The total CNDO/2 energy is: " ;
                totalEnergy = (nre + ee_alpha + ee_beta);
               cout << totalEnergy << " eV" << endl;
                P_alpha_final = P_alpha;
                P_beta_final = P_beta;
                nuclear_repulsion_energy = nre;
                electronic_energy = ee_alpha + ee_beta;
            }
            else {
                cout << "SCF has not converged" << endl;
            }
        }
        iteration += 1;
    }

}

// Function: Solve X Matrix 
mat Molecule::solveXMatrix(mat P_alpha, mat P_beta) {
    mat X = zeros<mat>(numBasisSets, numBasisSets);

    // add alpha and beta density together
    mat newMatrix = P_alpha + P_beta;
    // cout << newMatrix << endl;

    int mu = 0;
    for (int A = 0; A < numAtoms; A++) {
        string name_A = atomList[A].atomName;
        int numBasis_A = atomList[A].basisFunctionList.size();

        // loop over all the basis sets in atomA
        for (int i = 0; i < numBasis_A; i++) {
            int nu = 0;
            // loop over all the atoms in the molecule 
            for (int B = 0; B < numAtoms; B++) {
                string name_B = atomList[B].atomName;
                int numBasis_B = atomList[B].basisFunctionList.size();
                // loop over all the basis sets in atomB 
                for (int j = 0; j < numBasis_B; j++) {
                    if (mu != nu) {
                        int beta_A = atom_para[name_A](2);
                        int beta_B = atom_para[name_B](2);
                        X(mu, nu) = (beta_A + beta_B) * newMatrix(mu, nu);
                    }
                    nu += 1;
                }
            }
            mu += 1;
        }
    }
    return X;
}

// Function: Solve the Y Matrix 
mat Molecule::solveYMatrix(mat P_alpha, mat P_beta) {
    mat Y = zeros<mat>(numAtoms, numAtoms);

    // determine the total density 
    vec totalDensityVec = totalDensity(P_alpha, P_beta);

    int mu = 0;
    // loop over all the atoms
    for (int A = 0; A < numAtoms; A++) {
        for (int B = 0; B < numAtoms; B++) {
            Y(A, B) = totalDensityVec(A) * totalDensityVec(B)
                - (atomList[A].valence * totalDensityVec(B) + atomList[B].valence * totalDensityVec(A));
            int vu = 0;
            for (int i = 0; i < atomList[A].basisFunctionList.size(); i++) {
                for (int j = 0; j < atomList[B].basisFunctionList.size(); j++) {
                    Y(A, B) -= ((P_alpha(mu, vu) * P_alpha(mu, vu)) + (P_beta(mu, vu) * P_beta(mu, vu)));
                }
                vu += 1;
            }
        }
        mu += 1;
    }
    return Y;
}


// Function: Calculate the contracted overlap integral for two primitive Gaussian shells gradient 
vec Molecule::contractedOverlapIntegral_derivative(basisFunction basisFunction1, basisFunction basisFunction2) {
    double overlap_dx = 0;
    double overlap_dy = 0;
    double overlap_dz = 0;

    vec center = zeros<vec>(3);
    vec result = zeros<vec>(3);

    int dim = 0;

    // loop over all the exponent combinations for dimension 1 
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            overlap_dx += basisFunction1.coefficients(i) * basisFunction2.coefficients(j) * basisFunction1.normalization(i) * basisFunction2.normalization(j) *
                calc_derivative_vals(0, basisFunction1, basisFunction2) *
                overlap_at_1D(basisFunction1.exponents(i), basisFunction2.exponents(j), basisFunction1.center(1), basisFunction2.center(1), basisFunction1.quantumNumbers(i), basisFunction2.quantumNumbers(j)) *
                overlap_at_1D(basisFunction1.exponents(i), basisFunction2.exponents(j), basisFunction1.center(2), basisFunction2.center(2), basisFunction1.quantumNumbers(i), basisFunction2.quantumNumbers(j));
        }
    }

    // loop over all the exponent combinations for dimension 2 
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            overlap_dy += basisFunction1.coefficients(i) * basisFunction2.coefficients(j) * basisFunction1.normalization(i) * basisFunction2.normalization(j) *
                calc_derivative_vals(1, basisFunction1, basisFunction2);
            overlap_at_1D(basisFunction1.exponents(i), basisFunction2.exponents(j), basisFunction1.center(0), basisFunction2.center(0), basisFunction1.quantumNumbers(i), basisFunction2.quantumNumbers(j))*
                overlap_at_1D(basisFunction1.exponents(i), basisFunction2.exponents(j), basisFunction1.center(2), basisFunction2.center(2), basisFunction1.quantumNumbers(i), basisFunction2.quantumNumbers(j));
        }
    }

    // loop over all the exponent combinations for dimension 2 
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            overlap_dy += basisFunction1.coefficients(i) * basisFunction2.coefficients(j) * basisFunction1.normalization(i) * basisFunction2.normalization(j) *
                calc_derivative_vals(2, basisFunction1, basisFunction2);
            overlap_at_1D(basisFunction1.exponents(i), basisFunction2.exponents(j), basisFunction1.center(0), basisFunction2.center(0), basisFunction1.quantumNumbers(i), basisFunction2.quantumNumbers(j))*
                overlap_at_1D(basisFunction1.exponents(i), basisFunction2.exponents(j), basisFunction1.center(1), basisFunction2.center(1), basisFunction1.quantumNumbers(i), basisFunction2.quantumNumbers(j));
        }
    }

    result(0) = overlap_dx;
    result(1) = overlap_dy;
    result(2) = overlap_dz;

    return result;
}

// Function: Calculate the overlap integral for two primitive Gaussian shells gradient 
arma::field<arma::vec> Molecule::overlapMatrix_derivative() {
    arma::field<arma::vec> overlap(numBasisSets, numBasisSets);

    int mu = 0;
    for (int A = 0; A < numAtoms; A++) {
        int numBasis_A = atomList[A].basisFunctionList.size();
        for (int i = 0; i < numBasis_A; i++) {
            int nu = 0;
            for (int B = 0; B < numAtoms; B++) {
                int numBasis_B = atomList[B].basisFunctionList.size();
                for (int j = 0; j < numBasis_B; j++) {
                    overlap(mu, nu) = contractedOverlapIntegral_derivative(atomList[A].basisFunctionList[i], atomList[B].basisFunctionList[j]);
                    nu += 1;
                }
            }
            mu += 1;
        }
    }
    return overlap;
}

// Function: Calculate the derivative value for the overlap integral
double Molecule::calc_derivative_vals(int dim, basisFunction basisFunction1, basisFunction basisFunction2) {
    // get the lmn values 
    vec lmn1 = basisFunction1.quantumNumbers;
    vec lmn2 = basisFunction2.quantumNumbers;

    double result;

    if (lmn1(dim) == 0) {
        result = 2 * basisFunction1.exponents(dim) * overlap_at_1D(basisFunction1.exponents(dim), basisFunction2.exponents(dim), basisFunction1.center(dim), basisFunction2.center(dim), basisFunction1.quantumNumbers(dim) + 1, basisFunction2.quantumNumbers(dim));
    }

    if (lmn1(dim) == 1) {
        result = (-1.0) * overlap_at_1D(basisFunction1.exponents(dim), basisFunction2.exponents(dim), basisFunction1.center(dim), basisFunction2.center(dim), basisFunction1.quantumNumbers(dim) - 1, basisFunction2.quantumNumbers(dim)) +
            2 * basisFunction1.exponents(dim) * overlap_at_1D(basisFunction1.exponents(dim), basisFunction2.exponents(dim), basisFunction1.center(dim), basisFunction2.center(dim), basisFunction1.quantumNumbers(dim) + 1, basisFunction2.quantumNumbers(dim));
    }

    return result;
}

// Function: Calculate the 0 to 0 value for Gamma Matrix
vec Molecule::calc_zero_to_zero(basisFunction basisFunction1, basisFunction basisFunction2, double first_term, double second_term) {

    vec result = zeros<vec>(3);

    double dist = norm(basisFunction1.center - basisFunction2.center);
    double v_2 = 1.0 / (first_term + second_term);
    double srT = sqrt(v_2) * dist;
    double T = pow(srT, 2);
    double u = pow(M_PI * first_term, 1.5) * pow(M_PI * second_term, 1.5);

    if (dist != 0.0) {
        result = (u) * (1.0 / pow(dist, 2)) * ((2 * sqrt(v_2) * exp(-T)) / sqrt(M_PI) - erf(srT) / dist) * (basisFunction1.center - basisFunction2.center);
    }

    return result * Hartree;
}

// Function: Calculate the derivative of the elements in the gamma matrix 
vec Molecule::calcGammaDerivative(basisFunction basisFunction1, basisFunction basisFunction2) {
    vec coeff_a = zeros<vec>(3);
    vec coeff_b = zeros<vec>(3);

    vec alpha_a = basisFunction1.exponents;
    vec alpha_b = basisFunction2.exponents;

    // d prime 
    for (int i = 0; i < 3; i++) {
        coeff_a(i) = basisFunction1.coefficients(i) * basisFunction1.normalization(i);
        coeff_b(i) = basisFunction2.coefficients(i) * basisFunction2.normalization(i);
    }

    vec val = zeros<vec>(3);

    for (int k = 0; k < 3; k++) {
        for (int j = 0; j < 3; j++) {
            double first_term = 1.0 / (alpha_a(k) + alpha_a(j));

            for (int l = 0; l < 3; l++) {
                for (int m = 0; m < 3; m++) {
                    double second_term = 1.0 / (alpha_b(l) + alpha_b(m));
                    val += coeff_a(k) * coeff_a(j) * coeff_b(l) * coeff_b(m) *
                        calc_zero_to_zero(basisFunction1, basisFunction2, first_term, second_term);
                }
            }
        }
    }
    return val;
}

// Function: Calculate and loop through the elements for Gamma Derivative Matrix 
field<vec> Molecule::calcGammaDerivativeMatrix() {

    field<vec> gammaDerivative(numAtoms, numAtoms);

    for (int A = 0; A < numAtoms; A++) {
        for (int B = 0; B < numAtoms; B++) {
            gammaDerivative(A, B) = calcGammaDerivative(atomList[A].basisFunctionList[0], atomList[B].basisFunctionList[0]);
        }
    }
    return gammaDerivative;
}

// Function: Calculate the gradient repulsion energy 
mat Molecule::calcGradientRepulsionEnergy() {
    mat gradient = zeros<mat>(numAtoms, 3);

    for (int A = 0; A < numAtoms; A++) {
        for (int B = 0; B < numAtoms; B++) {
            if (A != B) {
                vec3 delta = atomList[A].coords - atomList[B].coords;
                double dist = norm(delta);

                const double epsilon = 1e-6;
                double denominator = pow(dist, 3) + epsilon;


                // Update the gradient matrix
                gradient.row(A) += (atomList[A].valence * atomList[B].valence) / denominator * trans(delta);

            }
        }
    }

    // transpose the gradient 
    mat final_result = zeros<mat>(numAtoms, 3);
    final_result = gradient.t();

    return -1 * (final_result * Hartree);
}

// Function: Calculate the final gradient value 
mat Molecule::calcFinalGradient() {
    mat gradient = zeros<mat>(numAtoms, 3);

    // calculate the X and Y matrices
    mat X = solveXMatrix(P_alpha_final, P_beta_final);
    mat Y = solveYMatrix(P_alpha_final, P_beta_final);

    // calculate the gamma derivative matrix
    field<vec> gammaDerivative = calcGammaDerivativeMatrix();

    // calculate the overlap derivative matrix
    field<vec> overlapDerivative = overlapMatrix_derivative();

    // calculate the gradient repulsion energy
    mat gradientRepulsion = calcGradientRepulsionEnergy();

    int mu = 0;

    for (int A = 0; A < numAtoms; A++) {
        vec sum = zeros<vec>(3);
        int numBasis_A = atomList[A].basisFunctionList.size();
        for (int i = 0; i < numBasis_A; i++) {
            int nu = 0;
            for (int B = 0; B < numAtoms; B++) {
                int numBasis_B = atomList[B].basisFunctionList.size();
                for (int j = 0; j < numBasis_B; j++) {
                    if (mu != nu) {
                        sum += (X(mu, nu) * overlapDerivative(mu, nu));
                    }
                    if (A != B) {
                        sum += (Y(A, B) * gammaDerivative(A, B));
                    }
                }
                nu += 1;
            }
        }
        mu += 1;

        // add the gradient repulsion energy
        sum += gradientRepulsion.col(A);
        gradient.row(A) = sum.t();
    }

    return gradient.t();
}

