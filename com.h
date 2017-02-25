#include <iostream>
#include <complex>
#include <Eigen/Dense>
#include "OP.h"
#include "JC_Parameter.h"
#include "Hamiltanian.h"
#include "Evolution.h"

using namespace std;
using namespace Eigen;
void com(const JC_Parameter&para)
{
        
        Hamiltanian Qubit(qubit, para);
        Hamiltanian Resonator(resonator, para);

        Hamiltanian Ham1(Qubit), Ham2;

        

        for(int i=0; i<para.LatticeSize(); ++i)
        {
                if(i%2==0)
                        Ham2.kron(Ham1, Resonator, para.gr());
                else
                        Ham1.kron(Ham2, Qubit, para.gl());
        }


        MatrixXd H1=Ham1.System().QMat()->at(para.ParticleNo());

        cout<<H1.rows()<<"x"<<H1.cols()<<endl;

        Ham1=Qubit;
        for(int i=0; i<para.LatticeSize(); ++i)
        {
                if(i%2==0)
                        Ham2.kron(Ham1, Resonator, para.gl());
                else
                        Ham1.kron(Ham2, Qubit, para.gr());
        }

        MatrixXd H2=Ham1.System().QMat()->at(para.ParticleNo());
        cout<<H2.rows()<<"x"<<H2.cols()<<endl;

        Evolution tevo(H1, H2, 1);

        VectorXcd a(tevo._tOP*tevo._eigenstate.eigenvectors().col(1));

        cout<<a.adjoint()*tevo._eigenstate.eigenvectors().col(1)<<endl;



}