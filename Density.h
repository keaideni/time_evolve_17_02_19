//an array to store the density

#ifndef DENSITY_H
#define DENSITY_H
#include "OP.h"
#include "JC_Parameter.h"
#include "Hamiltanian.h"

class Density
{
private:
        std::vector<MatrixXd> _array;
public:
        Density(){};
        ~Density(){};
        Density(const Density& den):_array(den._array){};
        Density(const JC_Parameter& para);

        const std::vector<MatrixXd>& array()const{return _array;};

        
};



#endif // DENSITY_H


Density::Density(const JC_Parameter& para)
{
        

        Hamiltanian Qubit(qubit, para), Resonator(resonator, para);

        OP temp1, temp2;

        for(int i=0; i<para.LatticeSize(); ++i)
        {

                if(i==0)
                {
                        temp1.time(Qubit.SysCDagR(), Qubit.SysCR());
                        for(int j=1; j<para.LatticeSize(); ++j)
                        {
                                if(j%2==1)temp2.kron(temp1, Resonator.SysEye());
                                else temp2.kron(temp1, Qubit.SysEye());

                                temp1=temp2;
                        }

                                

                        
                        

                }else
                {
                        temp1=Qubit.SysEye();
                        for(int j=1; j<para.LatticeSize(); ++j)
                        {
                                if(j%2==1)
                                {
                                        if(j==i)
                                                temp2.kron(temp1, Resonator.SysCDagR()*Resonator.SysCR());
                                        else temp2.kron(temp1, Resonator.SysEye());
                                }else 
                                {
                                        if(j==i)temp2.kron(temp1, Qubit.SysCDagR()*Qubit.SysCR());
                                        else temp2.kron(temp1, Qubit.SysEye());
                                }
                                temp1=temp2;
                        }

                                
                }

                _array.push_back(temp1.QMat().at(para.ParticleNo()));
        }

}