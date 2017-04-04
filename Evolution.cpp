#include "Evolution.h"
#include <iostream>

using namespace std;



Evolution::Evolution(const MatrixXd& H1, const MatrixXd& H2, const double& t):
_tepOP(H2)
{
        

        //cout<<eigenstate.eigenvalues().size()<<endl;
        //cout<<tepOP.eigenvalues().size()<<endl;

        int n=_tepOP.eigenvalues().size();
        VectorXcd timeOP(n);

        for(int i=0; i<n; ++i)timeOP[i]=std::complex<double>(cos(_tepOP.eigenvalues()(i)*t),-1*sin(_tepOP.eigenvalues()(i)*t));

        MatrixXcd E=timeOP.asDiagonal();
        //cout<<E.adjoint()*E<<endl;
        //cout<<"the real:"<<endl;
        //cout<<tepOP.eigenvalues()<<endl;
        //cout<<"the complex"<<endl;
        //cout<<timeOP<<endl;
        //cout<<"the diagonal:"<<endl;
        //cout<<E<<endl;
        
        //cout<<tepOP.eigenvectors()*tepOP.eigenvectors().inverse()<<endl;;
        _tOP=_tepOP.eigenvectors()*E*_tepOP.eigenvectors().adjoint();
        //cout<<_tOP<<endl;
        for(int i=0; i<n; ++i)timeOP[i]=std::complex<double>(1,0);
        _t0OP=(timeOP.asDiagonal());

        


//===========calculation of the initial state===============================
        DenseSymMatProd<double> op(H1);
        SymEigsSolver< double, SMALLEST_ALGE, DenseSymMatProd<double> > eigs(&op, 2, 6);
        eigs.init();
        eigs.compute();

        //if(eigs.info() == SUCCESSFUL)
        _eigenstate = eigs.eigenvectors().col(1);//change here to change the evolve state.

        
}