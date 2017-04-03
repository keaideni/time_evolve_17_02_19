#include <iostream>
#include <iostream>
#include <complex>
#include <Eigen/Dense>
#include "Mat.h"
#include "OP.h"
#include "JC_Parameter.h"
#include "Hamiltanian.h"
#include "Evolution.h"
#include "CalQ.h"
#include "mpi.h"

int OP::Max;

int main(int argc, char* argv[])
{
        ifstream inpara("./data/QNosave.txt");
        if(!inpara)
        {
                cerr<<"the file QNosave.txt doesn't exit!"<<endl;
        }

        JC_Parameter para(inpara);

        inpara.close();

        //cout<<"the Parameter is "<<endl;
        //para.show();


        


        OP::Max=para.ParticleNo();


        MPI_Status status;

        int myid, numprocess;

        int groupn(1000);
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocess);


        /*create a MPI_Datatype to express the complex<double>*/
        complex<double> car;
        const int nitems=2;
        int          blocklengths[2] = {1,1};
        MPI_Datatype types[2] = {MPI_LONG_DOUBLE, MPI_LONG_DOUBLE};
        MPI_Datatype mpi_car_type;
        MPI_Aint     offsets[2];

        offsets[0] = car.real();
        offsets[1] = car.imag();

        MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_car_type);
        MPI_Type_commit(&mpi_car_type);


        int everygroup(groupn/numprocess);//cout<<everygroup;

        Mat H1(para, 0.008, 0.022);
        Mat H2(para, 0.022, 0.008);

        

        //std::cout<<den.array()[1];

        if(myid==0)
        {

                Density den(para);

                ofstream qubitout("./result/qubit"), resonatorout("./result/resonator");
                ofstream entout("./result/entanglement");
                ofstream difout("./result/difference");
                ofstream energyout("./result/energy");


                Evolution tevo(H1.MatH(), H2.MatH(), 1);
                VectorXcd a=tevo.t0OP()*tevo.eigenstate();
                int Msize(tevo.tOP().cols());

                MPI_Bcast(&Msize, 1, MPI_INT, 0, MPI_COMM_WORLD);
                {
                        MatrixXcd temptop(tevo.tOP());
                        MPI_Bcast(&temptop(0,0), Msize*Msize, mpi_car_type, 0, MPI_COMM_WORLD);
                        
                }
                //=============================first time==================================

                outfile(den, a, tevo,  H2, difout, energyout, qubitout, resonatorout, entout, para);


                //==========================================================================
                for(int t=1; t<=everygroup; ++t)
                {
                        outfile(den, a, tevo,  H2, difout, energyout, qubitout, resonatorout, entout, t, para);
                }

                for(int id=1; id<numprocess; ++id)
                {
                        VectorXi trec(everygroup);
                        VectorXd energyrec(everygroup), entrec(everygroup);
                        MatrixXd difrec(everygroup,2), qubitdenrec(everygroup, para.ParticleNo());
                        MatrixXd resonatordenrec(everygroup, para.ParticleNo());

                        MPI_Recv(&trec(0), everygroup, MPI_INT, id, id, MPI_COMM_WORLD, &status);
                        MPI_Recv(&energyrec(0), everygroup, MPI_DOUBLE, id, id+numprocess,
                                        MPI_COMM_WORLD, &status);
                        MPI_Recv(&entrec(0), everygroup, MPI_DOUBLE, id, id+2*numprocess, 
                                        MPI_COMM_WORLD, &status);
                        MPI_Recv(&difrec(0,0), everygroup*2, MPI_DOUBLE, id, id+3*numprocess, 
                                        MPI_COMM_WORLD, &status);
                        MPI_Recv(&qubitdenrec(0,0), everygroup*para.ParticleNo(), MPI_DOUBLE, id,
                                        id+4*numprocess, MPI_COMM_WORLD, &status);
                        MPI_Recv(&resonatordenrec(0,0), everygroup*para.ParticleNo(), MPI_DOUBLE, id,
                                        id+5*numprocess, MPI_COMM_WORLD, &status);

                        outfile(trec, energyrec, entrec, difrec, qubitdenrec, resonatordenrec, difout, energyout,
                                qubitout, resonatorout, entout, para, everygroup);

                }

                energyout.close();entout.close();difout.close();qubitout.close();resonatorout.close();

        }else
        {

                Density den(para);

                int Msize;
                MPI_Bcast(&Msize, 1, MPI_INT, 0, MPI_COMM_WORLD);

                MatrixXcd tOP(Msize, Msize);
                MPI_Bcast(&tOP(0,0), Msize*Msize, mpi_car_type, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&tOP.imag()(0,0), Msize*Msize, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);


                //std::cout<<tOP<<endl;

                VectorXi t(everygroup);

                VectorXd energysend(everygroup), entsend(everygroup);
                MatrixXd difsend(everygroup,2), qubitdensend(everygroup, para.ParticleNo());
                MatrixXd resonatordensend(everygroup, para.ParticleNo());

                Evolution tevo(H1.MatH(), H2.MatH(), myid*everygroup);
                
                calculate(para, t, energysend, entsend, difsend, qubitdensend, resonatordensend, tOP, tevo,
                                H2, den, myid, everygroup);

                MPI_Send(&t(0), everygroup, MPI_INT, 0, myid, MPI_COMM_WORLD);
                MPI_Send(&energysend(0), everygroup, MPI_DOUBLE, 0, myid+numprocess, MPI_COMM_WORLD);
                MPI_Send(&entsend(0), everygroup, MPI_DOUBLE, 0, myid+2*numprocess, MPI_COMM_WORLD);
                MPI_Send(&difsend(0,0), everygroup*2, MPI_DOUBLE, 0, myid+3*numprocess, MPI_COMM_WORLD);
                MPI_Send(&qubitdensend(0,0), everygroup*para.ParticleNo(), MPI_DOUBLE, 0,
                         myid+4*numprocess, MPI_COMM_WORLD);
                MPI_Send(&resonatordensend(0,0), everygroup*para.ParticleNo(), MPI_DOUBLE, 0,
                                        myid+5*numprocess, MPI_COMM_WORLD);



        }


        MPI_Finalize();
        
        
}