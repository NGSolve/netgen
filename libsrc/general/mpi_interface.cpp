/**************************************************************************/
/* File:   mpi_interface.cpp                                              */
/* Author: Joachim Schoeberl                                              */
/* Date:   04. Apr. 97                                                    */
/**************************************************************************/

#include <mystdlib.h>
#include <myadt.hpp>


namespace netgen
{


#ifdef PARALLEL
  
  void MyMPI_SendCmd (const char * cmd)
  {
    int ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    if(ntasks==1)
      return;
    
    for (int dest = 1; dest < ntasks; dest++)
      MPI_Send( cmd, (strlen(cmd)+1), MPI_CHAR, dest, MPI_TAG_CMD, MPI_COMM_WORLD);
  }

  string MyMPI_RecvCmd ()
  {
    MPI_Status status;
    int flag;
    int size_of_msg = -1;

    MPI_Probe(0, MPI_TAG_CMD, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_CHAR, &size_of_msg);

    //char* buf = (char*)malloc(size_of_msg*sizeof(char));
    char buf[100000]; //1MB should be enough...

    MPI_Recv( &buf, size_of_msg, MPI_CHAR, 0, MPI_TAG_CMD, MPI_COMM_WORLD, &status);
    
    return string(buf);
  }

#endif


}

