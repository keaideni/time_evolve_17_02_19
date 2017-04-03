#include "OP.h"
#include "JC_Parameter.h"
#include "Hamiltanian.h"
#include "Density.h"
#include <fstream>

using namespace std;

void CalParNo(const Density& den, const VectorXcd& wave, VectorXd &qubitout, VectorXd& resonatorout, const Parameter& para);
void CalParNo(const Density& den, const VectorXcd& wave, VectorXd &qubitout, VectorXd& resonatorout, const Parameter& para)
{
        
 
        for(int i=0; i<para.LatticeSize(); ++i)
        {
                //cout<<i<<endl;

                if(i%2==0)
                {       
                        //std::cout<<(wave.adjoint()*den.array()[i]*wave).real()<<endl;
                        
                        qubitout(i/2)=(wave.adjoint()*den.array()[i]*wave)(0,0).real();
                }
                else
                        resonatorout(i/2)=(wave.adjoint()*den.array()[i]*wave)(0,0).real();

                       
                
        }

}



double calEnt(const JC_Parameter& para, const VectorXcd& wave);
double calEnt(const JC_Parameter& para, const VectorXcd& wave)
{
        Hamiltanian Qubit(qubit, para), Resonator(resonator, para);

        OP temp1(Qubit.SysEye()), temp2;
        for(int i=1; i<para.ParticleNo(); ++i)
        {
                if(i%2==1)temp2.kron(temp1, Resonator.SysEye());
                else temp2.kron(temp1, Qubit.SysEye());

                temp1=temp2;
        }
        OP saveOP(temp1);

        if(para.ParticleNo()%2==0)temp1=Qubit.SysEye();
        else temp1=Resonator.SysEye();

        for(int i=para.ParticleNo()+1; i<para.LatticeSize(); ++i)
        {
                if(i%2==1)temp2.kron(temp1, Resonator.SysEye());
                else temp2.kron(temp1, Qubit.SysEye());

                temp1=temp2;
        }

        OP waveOP;
        waveOP.Waveinitial(saveOP, temp1, para.ParticleNo());
        int i(0);
        double ent(0);
        for(auto it=waveOP.QMat().rbegin(); it!= waveOP.QMat().rend(); ++it)
        {
                MatrixXcd temp(it->second.rows(),it->second.cols());

                for(int rown=0; rown<it->second.rows(); ++rown)
                {
                        for(int coln=0; coln<it->second.cols(); ++coln)
                        {
                                temp(rown, coln)=wave(i++);
                                if(i>wave.size())
                                {
                                        cerr<<"the wave size is wrong!"<<endl;
                                        exit(1);
                                }
                        }
                }

                JacobiSVD<MatrixXcd> svd(temp, ComputeThinU | ComputeThinV);

                for(int j=0; j<svd.singularValues().size(); ++j)
                {
                        ent-=pow(svd.singularValues()(j),2)*log(pow(svd.singularValues()(j),2));
                }
        }



        return ent;

}

//the t=0 one, which don't need top to multiply.
void outfile(const Density& den, VectorXcd& a, const Evolution& tevo, const Mat& H2,
                 ofstream& difout, ofstream& energyout, ofstream& qubitout, 
                 ofstream& resonatorout, ofstream& entout, const JC_Parameter& para);
void outfile(const Density& den, VectorXcd& a, const Evolution& tevo, const Mat& H2,
                 ofstream& difout, ofstream& energyout, ofstream& qubitout, 
                 ofstream& resonatorout, ofstream& entout, const JC_Parameter& para)
{
        

        //=====================for the difference between waves====================
                
        difout<<0<<"\t"<<abs((a.adjoint()*tevo.eigenstate())(0,0))
        <<"\t"<<abs((a.adjoint()*tevo.tepOP().eigenvectors().col(0))(0,0))<<endl;
        //=========================================================================

        //===========energy===================
        energyout<<0<<"\t"<<(a.adjoint()*H2.MatH()*a).real()<<endl;

        //cout<<tevo._tOP.adjoint()*H2*tevo._tOP<<endl<<H2<<endl;
        VectorXd qubitden(para.LatticeSize()/2), resonatorden(para.LatticeSize()/2);

        CalParNo(den, a, qubitden, resonatorden, para);
        qubitout<<0<<"\t";resonatorout<<0<<"\t";
        for(int i=0; i<para.ParticleNo(); ++i)
        {
                qubitout<<qubitden(i)<<"\t";
                resonatorout<<resonatorden(i)<<"\t";
        }
        qubitout<<endl;resonatorout<<endl;

        entout<<0<<"\t"<<calEnt(para, a)<<endl;
}





void outfile(const Density& den, VectorXcd& a, const Evolution& tevo, const Mat& H2,
                 ofstream& difout, ofstream& energyout, ofstream& qubitout, 
                 ofstream& resonatorout, ofstream& entout, const int& t, const JC_Parameter& para);
void outfile(const Density& den, VectorXcd& a, const Evolution& tevo, const Mat& H2,
                 ofstream& difout, ofstream& energyout, ofstream& qubitout, 
                 ofstream& resonatorout, ofstream& entout, const int& t, const JC_Parameter& para)
{
        VectorXcd b(tevo.tOP()*a);
        a=b;

        //=====================for the difference between waves====================
                
        difout<<t<<"\t"<<abs((a.adjoint()*tevo.eigenstate())(0,0))
        <<"\t"<<abs((a.adjoint()*tevo.tepOP().eigenvectors().col(0))(0,0))<<endl;
        //=========================================================================

        //===========energy===================
        energyout<<t<<"\t"<<(a.adjoint()*H2.MatH()*a).real()<<endl;

        //cout<<tevo._tOP.adjoint()*H2*tevo._tOP<<endl<<H2<<endl;
        VectorXd qubitden(para.LatticeSize()/2), resonatorden(para.LatticeSize()/2);

        CalParNo(den, a, qubitden, resonatorden, para);
        qubitout<<t<<"\t";resonatorout<<t<<"\t";
        for(int i=0; i<para.ParticleNo(); ++i)
        {
                qubitout<<qubitden(i)<<"\t";
                resonatorout<<resonatorden(i)<<"\t";
        }
        qubitout<<endl;resonatorout<<endl;

        entout<<t<<"\t"<<calEnt(para, a)<<endl;
}


void outfile(const VectorXi& trec, const VectorXd& energyrec, const VectorXd& entrec,
                const MatrixXd& difrec, const MatrixXd& qubitdenrec, const MatrixXd& resonatordenrec,
                ofstream& difout, ofstream& energyout, ofstream& qubitout, 
                 ofstream& resonatorout, ofstream& entout, const Parameter& para, const int& everygroup);
void outfile(const VectorXi& trec, const VectorXd& energyrec, const VectorXd& entrec,
                const MatrixXd& difrec, const MatrixXd& qubitdenrec, const MatrixXd& resonatordenrec,
                ofstream& difout, ofstream& energyout, ofstream& qubitout, 
                 ofstream& resonatorout, ofstream& entout, const Parameter& para, const int& everygroup)
{
        for(int i=0; i<everygroup; ++i)
        {
                energyout<<trec(i)<<"\t"<<energyrec(i)<<endl;
                
                entout<<trec(i)<<"\t"<<entrec(i)<<endl;

                difout<<trec(i)<<"\t";
                for(int j=0; j<2; ++j)difout<<difrec(i,j)<<"\t";
                difout<<endl;

                qubitout<<trec(i)<<"\t";
                for(int j=0; j<para.LatticeSize()/2; ++j)qubitout<<qubitdenrec(i,j)<<"\t";
                qubitout<<endl;

                resonatorout<<trec(i)<<"\t";
                for(int j=0; j<para.LatticeSize()/2; ++j)resonatorout<<qubitdenrec(i,j)<<"\t";
                resonatorout<<endl;
        }
}


void calculate(const JC_Parameter& para, VectorXi& t,VectorXd& energysend,VectorXd& entsend, MatrixXd& difsend,
                MatrixXd& qubitdensend, MatrixXd& resonatordensend, const MatrixXcd& tOP,
                const Evolution& tevo, const Mat& H2, const Density& den, const int& myid, const int&everygroup);
void calculate(const JC_Parameter& para, VectorXi& t,VectorXd& energysend,VectorXd& entsend, MatrixXd& difsend,
                MatrixXd& qubitdensend, MatrixXd& resonatordensend, const MatrixXcd& tOP,
                const Evolution& tevo, const Mat& H2, const Density& den, const int& myid, const int&everygroup)
{
        VectorXcd a(tevo.tOP()*tevo.eigenstate());
        for(int i=0; i<everygroup; ++i)
        {
                VectorXcd b(tOP*a);
                a=b;

                t(i)=myid*everygroup+1+i;
                
                difsend(i,0)=abs((a.adjoint()*tevo.eigenstate())(0,0));
                difsend(i,1)=abs((a.adjoint()*tevo.tepOP().eigenvectors().col(0))(0,0));
                //=========================================================================

                //===========energy===================
                energysend(i)=(a.adjoint()*H2.MatH()*a)(0,0).real();

                //cout<<tevo._tOP.adjoint()*H2*tevo._tOP<<endl<<H2<<endl;
                VectorXd qubitden(qubitdensend.cols()), resonatorden(resonatordensend.cols());

                CalParNo(den, a, qubitden, resonatorden, para);
                
                qubitdensend.row(i)=qubitden.transpose();
                resonatordensend.row(i)=resonatorden.transpose();

                entsend(i)=calEnt(para, a);
        }
}