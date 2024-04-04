#include "molecule.h"
#include "utils.h"
#include <iostream>

using namespace std;

int main() {
    string filename = "H2.txt";
    Molecule mol(filename);

    double convergence = 1e-6;
    int maxIterations = 100;
    int situation = 0;
    mol.runSCF(convergence, maxIterations,situation);

    cout << "Suv_RA (A is the center of u)" << endl;
    print_field(mol.overlapMatrix_derivative());

    cout << "Gamma_AB_RA;" << endl;
    cout << mol.calcGammaDerivativeMatrix() << endl;

    cout << "gradient (Nuclear part)" << endl;
    cout << mol.calcGradientRepulsionEnergy() << endl;

    cout << "gradient (Electron part)" << endl;
    cout << mol.calcFinalGradient() - (mol.calcGradientRepulsionEnergy()) << endl;

    cout << "gradient" << endl;
    cout << mol.calcFinalGradient() << endl;

    printf("\nIf I use finite difference, change the x coordinate of second atom:\n");
    
    cout << "The molecule in the input file(Situation I)" << filename << endl;
;
    mol.runSCF(convergence, maxIterations,1);
    //Number is multiplied with 1e-6 to change the x coordinates, if you wanna change the y coordinates,
    //just to change the atomlist[0 to 1].coord.
    //    atomList[1].coords(0) = atomList[1].coords(0) + situation * 1e-6;

    cout << "The molecule in the input file(Situation II)" << filename << endl;
    mol.runSCF(convergence, maxIterations,-1);

    cout << "The molecule in the input file(Situation III)" << filename<<endl;
    mol.runSCF(convergence, maxIterations,0);
    cout << "The molecule has the gradiant like this:" << endl;
    cout << mol.calcFinalGradient() << endl;

    return 0;
}