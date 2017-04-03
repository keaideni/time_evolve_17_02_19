#include <iostream>
#include <complex>
#include <Eigen/Dense>
#include "Mat.h"
#include "OP.h"
#include "JC_Parameter.h"
#include "Hamiltanian.h"
#include "Evolution.h"
#include "CalQ.h"

using namespace std;
using namespace Eigen;
void com(const JC_Parameter&para)
{
        
        Mat H1(para, para.gl(), para.gr());
        Mat H2(para, para.gr(), para.gl());




        //cout<<"first"<<endl;
        //cout<<H2.rows()<<"x"<<H2.cols()<<endl;

        ofstream qubitout("./result/qubit"), resonatorout("./result/resonator");
        ofstream entout("./result/entanglement");
        ofstream difout("./result/difference");
        ofstream energyout("./result/energy");


        Evolution tevo(H1.MatH(), H2.MatH(), 1);



        VectorXcd a=tevo._t0OP*tevo._eigenstate;

        VectorXcd b;

//=====================================first time======================================================



        //=====================for the difference between waves====================
                
        difout<<abs((a.adjoint()*tevo._eigenstate)(0,0))
        <<"\t"<<abs((a.adjoint()*tevo._tepOP.eigenvectors().col(0))(0,0))<<endl;
        //=========================================================================

        //===========energy===================
        energyout<<(a.adjoint()*H2.MatH()*a).real()<<endl;

        //cout<<tevo._tOP.adjoint()*H2*tevo._tOP<<endl<<H2<<endl;

        Hamiltanian Qubit(qubit, para);
        Hamiltanian Resonator(resonator, para);


        CalParNo(Qubit, Resonator, para, a, qubitout, resonatorout, para);

        calEnt(Qubit, Resonator, para, a, entout);
//==================================================================================================

        for(double t=1; t<1000; t+=1)
        {
                
                //cout<<tevo._tOP<<endl;

                 b=(tevo._tOP*a);

                 a=b;
//=====================for the difference between waves====================
                
                difout<<abs((a.adjoint()*tevo._eigenstate)(0,0))
                <<"\t"<<abs((a.adjoint()*tevo._tepOP.eigenvectors().col(0))(0,0))<<endl;
//=========================================================================

                //===========energy===================
                energyout<<abs((a.adjoint()*H2.MatH()*a)(0,0))<<endl;


                CalParNo(Qubit, Resonator, para, a, qubitout, resonatorout, para);

                calEnt(Qubit, Resonator, para, a, entout);

        }

        qubitout.close();resonatorout.close();
        entout.close();
        difout.close();
        energyout.close();


        //cout<<a.adjoint()*tevo._eigenstate.eigenvectors().col(1)<<endl;



}