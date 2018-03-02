''' File      : mpiVina.py
*   Author    : Jeremy Hofer and Abhay Aradhya
*   Purpose   : Python MPI based parallel version of Autodock Vina.
'''
from mpi4py import MPI
from collections import deque

#define MASTER                  0

#define COMPUTE_TAG             11
#define TERMINATE_TAG           22
#define WORK_REQ_TAG            33

#define MAX_LIGAND_NAME_LENGTH  25
#define LIGAND_FILE_NAME        "ligandlist"

MASTER = 0

COMPUTE_TAG = 11
TERMINATE_TAG = 22
WORK_REQ_TAG = 33

def main():
    numProcs = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    totalLigands = 0
    queue = deque(None)

    #Only master processor will read the ligandlist file and will make the work pool.
    if rank == MASTER:
        print "Master processor : Reading ligandlist file and creating work pool.\n"
        startTime = MPI.Wtime(); #start timer.
        #insert opening config files and reading in lists of receptors, ligands, etc.
        '''
        if (NULL == ligandListFile)
        {
            printf("Couldn't open file %s for reading.\n", LIGAND_FILE_NAME);
            MPI_Abort(MPI_COMM_WORLD, 911); //Terminates MPI execution environment with error code 911.
            return 0;
        }
        '''

        #append items to the queue

    if rank == MASTER:
        mpiVinaManager(numProcs);   #Master processor will play the role of mpiVINA manager.
    else:
        mpiVinaWorker(rank);    #All other processors will play the role of mpiVINA worker.

    if rank == MASTER:
        endTime = MPI.Wtime(); #end timer.
        print "\n\n..........................................\n" 
        print "   Number of workers       = {0} \n".format(numProcs - 1)
        print "   Number of Ligands        = {0} \n".format(totalLigands)
        print "   Total time required     = {0} seconds.\n".format(endTime - startTime)
        print "..........................................\n\n"

def mpiVinaManager(numProcs):
    mStatus = MPI.Status()
    while len(queue) > 0:
        comm.recv(source=MPI.ANY_SOURCE, tag=WORK_REQ_TAG, status=mStatus)
        comm.send(buf=queue.popleft(), dest=mStatus.Get_source(), tag=COMPUTE_TAG)

    for i in range(numProcs):
        comm.recv(source=MPI.ANY_SOURCE, tag=WORK_REQ_TAG, status=mStatus)
        comm.send(buf=None, dest=mStatus.Get_source(), tag=TERMINATE_TAG)

def mpiVinaWorker(workerID):
    print "Worker %d has started.\n".format(workerID)
    wStatus = MPI.Status()

    comm.send(None, dest=0, tag=WORK_REQ_TAG)
    ligandName = comm.recv(source=0, tag=MPI.ANY_TAG, status=wStatus)

    while wStatus.Get_tag == COMPUTE_TAG:
        print "Worker = {0} : ligand {1} is processing...\n".format(workerID, ligandName)
        #insert vina command
        comm.send(buf=None, dest=0, tag=WORK_REQ_TAG)
        ligandName = comm.recv(source=0, tag=MPI.ANY_TAG, status=wStatus)

    if wStatus.Get_tag == TERMINATE_TAG:
        print "Worker {0} has terminated.\n".format(workerID)
    else:
        print "Worker {0} has received invalid Tag\n".format(workerID)

if __name__ == '__main__':
    main()