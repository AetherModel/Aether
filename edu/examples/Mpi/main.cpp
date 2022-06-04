

#include <iostream>
#include <vector>
#include "mpi.h"
#include <unistd.h>

int main() {

  MPI_Comm aether_comm;
  int nProcs;
  int iProc;
  int iTag;

  MPI_Init(NULL, NULL);
  aether_comm = MPI_COMM_WORLD;

  // Get the number of processes
  MPI_Comm_size(aether_comm, &nProcs);

  // Get the rank of the process
  MPI_Comm_rank(aether_comm, &iProc);
  
  //std::cout << "iProc " << iProc << " of " << nProcs << "\n";

  MPI_Request *requests = static_cast<MPI_Request*>(malloc(sizeof(MPI_Request)*nProcs));
  
  int nVars = 300;
  std::vector<std::vector<float>> all_data;
  std::vector<float> tmp;
  float value;
  int iP, iV, iReceiver, iSender;
  
  int iSize = nVars * sizeof(float);
  std::vector<float*> buffers;
  std::vector<float*> rbuffers;
  
  // send all :

  for (int iP = 0; iP < nProcs; iP++) {
    buffers.push_back(static_cast<float*>(malloc(iSize)));
    if (iProc == 0)
      std::cout << "in iP : " << iP << " ";
    for (int iV = 0; iV < nVars; iV++) {
      value = iP * 100 + iV * 2 + 1;
      buffers[iP][iV] = value;
      if (iProc == 0)
	std::cout << value << " ";
    }
    if (iProc == 0)
      std::cout << "\n";
    iTag = iProc * 10 + iP;
    iReceiver = iP;
    MPI_Isend(buffers[iP], iSize, MPI_BYTE,
              iReceiver, iTag, aether_comm, &requests[iP]);
  }

  // once the processor sends everything, then it can recieve everything...

  for (int iP = 0; iP < nProcs; iP++) {
    iTag = iP * 10 + iProc;
    iSender = iP;
    rbuffers.push_back(static_cast<float*>(malloc(iSize)));
    MPI_Recv(rbuffers[iP], iSize, MPI_BYTE, iSender, iTag, aether_comm,
             MPI_STATUS_IGNORE);
  }
  
  for (int iP = 0; iP < nProcs; iP++) {
    MPI_Wait(&requests[iP], MPI_STATUS_IGNORE);
  }

  MPI_Barrier(aether_comm);
  sleep(0.1);
  for (int iP = 0; iP < nProcs; iP++) {
    if (iProc == iP) {
      std::cout << "iProc, iP : " << iProc << " " << iP << " ";
      for (int iV = 0; iV < nVars; iV++) {
	std::cout << rbuffers[iP][iV] << " ";
      }
      std::cout << "\n";
    }
    MPI_Barrier(aether_comm);
    sleep(0.1);
  }
  
  MPI_Barrier(aether_comm);
  sleep(1.0);
  
  MPI_Finalize();

  sleep(0.25);
  
  return 0;
}
