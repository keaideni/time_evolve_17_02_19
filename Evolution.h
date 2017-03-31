//here this class provides the eigenstates of initial Hamiltanian and the time evolution operators.


#ifndef EVOLUTION_H
#define EVOLUTION_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <SymEigsSolver.h>  // Also includes <MatOp/DenseSymMatProd.h>
#include "Hamiltanian.h"

using namespace Eigen;
using namespace Spectra;

class Evolution
{
private:
        
public:
        Evolution(){};
        ~Evolution(){};

        VectorXd _eigenstate;
        SelfAdjointEigenSolver<MatrixXd> _tepOP;
        MatrixXcd _tOP;

        MatrixXcd _t0OP;//used for the initial wave, make it complex.




        Evolution(const MatrixXd& H1, const MatrixXd& H2, const double& t);
};











#endif // EVOLUTION_H
