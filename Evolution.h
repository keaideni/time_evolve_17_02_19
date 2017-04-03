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
        VectorXd _eigenstate;
        SelfAdjointEigenSolver<MatrixXd> _tepOP;
        MatrixXcd _tOP;

        MatrixXcd _t0OP;//used for the initial wave, make it complex.
public:
        Evolution(){};
        ~Evolution(){};

        const VectorXd& eigenstate() const {return _eigenstate;};
        const SelfAdjointEigenSolver<MatrixXd>& tepOP()const {return _tepOP;};
        const MatrixXcd& tOP() const {return _tOP;};
        const MatrixXcd& t0OP() const {return _t0OP;};




        Evolution(const MatrixXd& H1, const MatrixXd& H2, const double& t);
};











#endif // EVOLUTION_H
