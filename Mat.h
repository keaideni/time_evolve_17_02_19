//this one here is for the calculation of Hamiltanian matrix. Just calculate the matrix, we don't care about the
//eigenstates and eigenvalues here.




#ifndef MAT_H
#define MAT_H
#include <Eigen/Core>
#include <SymEigsSolver.h>  // Also includes <MatOp/DenseSymMatProd.h>
#include <iostream>
#include "Hamiltanian.h"

using namespace Eigen;
using namespace Spectra;
class Mat
{
private:
        MatrixXd _MatH;
        //MatrixXd _state;
        //VectorXd _Eigenvalues;
public:
        Mat(){};
        ~Mat(){};
        MatrixXd MatH()const{return _MatH;};
        //VectorXd state()const{return _state;};
        //VectorXd Eigenvalues()const{return _Eigenvalues;};



        Mat(const JC_Parameter& para, const double& gl, const double& gr)
        {
                Hamiltanian Qubit(qubit, para);
                Hamiltanian Resonator(resonator, para);
                //Qubit.show();Resonator.show();

                Hamiltanian Ham1(Qubit), Ham2;
                //Ham1.show();
                

                for(int i=1; i<para.LatticeSize(); ++i)
                {
                        //std::cout<<gl<<std::endl;
                        if(i%2==0)
                                Ham2.kron(Ham1, Qubit, gl);
                        else
                        {
                                Ham2.kron(Ham1, Resonator, gr);
                                //Ham2.show();
                        }
                        Ham1=Ham2;

                }
                Ham1.final(gl);


                //int coln(Ham1.System().QMat()->at(para.ParticleNo()).cols());
                //int rown(Ham1.System().QMat()->at(para.ParticleNo()).rows());

                //_MatH.resize(rown, coln);
                //Ham1.System().show();
                _MatH=Ham1.System().QMat()->at(para.ParticleNo());

                /*DenseSymMatProd<double> opmat(Ham1.System().QMat()->at(para.ParticleNo()));

                SymEigsSolver<double, SMALLEST_MAGN, DenseSymMatProd<double>>
                 eigs(&opmat, 5, 10);

                eigs.init();
                eigs.compute();

                _state=eigs.eigenvectors();
                _Eigenvalues=eigs.eigenvalues();*/


                /*SelfAdjointEigenSolver<MatrixXd> es(Ham1.System().QMat()->at(para.ParticleNo()));
                _state=es.eigenvectors();
                _Eigenvalues=es.eigenvalues();*/
        };
};








#endif