#include "Conductivity.hpp"

Conductivity::Conductivity(LSMSSystemParameters &lsms, LSMSCommunication &comm, LocalTypeInfo &local) {
    efermi = Complex(lsms.chempot,0.0001);
    cm.resize(local.atom.size(), 2);
    
    double timeSingleScatterers=MPI_Wtime();
    local.tmatStore=0.0;
    local.JxStore=0.0;
    local.JyStore=0.0;
    local.JzStore=0.0;
    expectTmatCommunication(comm,local);
    expectJxCommunication(comm,local);
    expectJyCommunication(comm,local);
    expectJzCommunication(comm,local); 
    
    if(lsms.global.iprint>=1) printf("calculate single scatterer solutions.\n");
    
    for (int i=0;i<local.atom.size();i++){
       for(int is=0;is<=0;is++){
          cm(i,is).init(lsms, local, local.atom[i], i, efermi, is);
       }
    }
    std::cout << std::endl;
    std::cout << "Atom 1: Current Mat" << std::endl;
    std::cout << std::endl;
    for (int i=0;i<16;i++){
       for (int j=0;j<16;j++){
          std::cout << local.JxStore(j+i*16,0) << "   ";
       }
       std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Atom 2: Current Mat" << std::endl;
    std::cout << std::endl;
    for (int i=0;i<16;i++){
       for (int j=0;j<16;j++){
          std::cout << local.JxStore(j+i*16,1) << "   ";
       }
       std::cout << std::endl;
    }
    std::cout << std::endl;

    if(lsms.global.iprint>=2) printf("About to send t matrices\n");
    sendTmats(comm,local);
    if(lsms.global.iprint>=2) printf("About to finalize t matrices communication\n");
    finalizeTmatCommunication(comm);
    if(lsms.global.iprint>=2) printf("Recieved all t matricies\n");

    if(lsms.global.iprint>=2) printf("About to send J matrices\n");
    sendJx(comm,local);
    sendJy(comm,local);
    sendJz(comm,local);
    if(lsms.global.iprint>=2) printf("About to finalize J matrices communication\n");
    finalizeTmatCommunication(comm);
    finalizeJxCommunication(comm);
    finalizeJyCommunication(comm);
    finalizeJzCommunication(comm);
    if(lsms.global.iprint>=2) printf("Recieved all J matricies\n");
    timeSingleScatterers=MPI_Wtime()-timeSingleScatterers;
    if(lsms.global.iprint>=0) printf("timeSingleScatteres = %lf sec\n",timeSingleScatterers);

}
