''' File      : mpiDOCK.py
*   Author    : Jeremy Hofer and Abhay Aradhya
*   Purpose   : Python MPI based parallel master/slave manager of Autodock Vina and various Vina rescoring functions
'''
from mpi4py import MPI
from collections import deque
from glob import glob

MASTER = 0

COMPUTE_TAG = 11
TERMINATE_TAG = 22
WORK_REQ_TAG = 33

config_file = "mpidock.config"

def abort_mpi(error_message):
    print error_message
    print "Terminating mpiDOCK"
    MPI.COMM_WORLD.Abort(911)
    exit()

def main():
    numProcs = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    configuration = {}

    if numProcs < 2:
        abort_mpi("Not enough processors! You must use at least 2 processors.")

    #Only master processor will read the ligandlist file and will make the work pool.
    if rank == MASTER:
        print "Master processor initializing mpiDOCK.".format(MASTER)
        print "Reading configuration file..."

        configuration = {"vina": None, "vina_config": None, "run_vina": False, "receptors_dir": None, "ligands_dir": None, "job_name": None,
                         "ligand_db_name": None, "output_dir": None, "clean_temp_files": None}
        
        #try to open and read in configuration file. abort MPI if errors occur
        try:
            file = open(config_file, "r")
        except (FileNotFoundError, IOError) as e:
            abort_mpi("Cannot find or open the configuration file. Please make sure mpidock.config is in the same directory as mpiDOCK.py.")
        else:
            lines = [x.strip().split("=") for x in file.readlines() if x[0] is not "#"]
            if len(lines) is 0:
                abort_mpi("Error reading configuration file.")

            for line in lines:
                if line[0] in configuration:
                    if line[0] in ["run_vina", "clean_temp_files"]:
                        if line[1] == "yes":
                            line[1] = True
                        else:
                            line[1] = False

                    configuration[line[0]] = line[1]
                else:
                    abort_mpi("Invalid key/variable {0} in configuration file. Please review valid keys/variables and try again.")

        #print configuration information
        print "Configuration read OK. Configured with the following parameters:"

        for key, value in configuration.iteritems():
            print "{0} : {1}".format(key, value)

        #initialize job queue
        queue = deque([])

        #append items to the queue
        for i in range(15):
            queue.append("test {0}".format(i))

        totalLigands = len(queue)

        print "---- STARTING mpiDOCK ----"

    MPI.COMM_WORLD.Barrier()
    MPI.COMM_WORLD.bcast(configuration, 0)

    if rank == MASTER:
        startTime = MPI.Wtime(); #start timer.
        mpiVinaManager(numProcs, queue);   #Master processor will play the role of mpiVINA manager.
    else:
        mpiVinaWorker(rank);    #All other processors will play the role of mpiVINA worker.

    MPI.COMM_WORLD.Barrier()

    if rank == MASTER:
        endTime = MPI.Wtime(); #end timer.
        print "\n.........................................." 
        print "   Number of workers       = {0}".format(numProcs - 1)
        print "   Number of Ligands        = {0}".format(totalLigands)
        print "   Total time required     = {0} seconds.".format(endTime - startTime)
        print ".........................................."

    MPI.Finalize()

def mpiVinaManager(numProcs, queue):
    print "Config in manager is {1}".format(workerID, configuration)
    mStatus = MPI.Status()
    while len(queue) > 0:
        MPI.COMM_WORLD.recv(source=MPI.ANY_SOURCE, tag=WORK_REQ_TAG, status=mStatus)
        print "Worker {0} requesting work\n".format(mStatus.Get_source())
        MPI.COMM_WORLD.send(queue.popleft(), dest=mStatus.Get_source(), tag=COMPUTE_TAG)

    for i in range(numProcs - 1):
        MPI.COMM_WORLD.recv(source=MPI.ANY_SOURCE, tag=WORK_REQ_TAG, status=mStatus)
        print "Sending termination to worker {0}\n".format(mStatus.Get_source())
        MPI.COMM_WORLD.send(None, dest=mStatus.Get_source(), tag=TERMINATE_TAG)

def mpiVinaWorker(workerID):
    print "Worker {0} has started.\n".format(workerID)
    print "Config in worker {0} is {1}".format(workerID, configuration)
    wStatus = MPI.Status()

    MPI.COMM_WORLD.send(None, dest=0, tag=WORK_REQ_TAG)
    ligandName = MPI.COMM_WORLD.recv(source=0, tag=MPI.ANY_TAG, status=wStatus)

    while wStatus.Get_tag() == COMPUTE_TAG:
        print "Worker = {0} : ligand {1} is processing...\n".format(workerID, ligandName)
        #insert vina command
        MPI.COMM_WORLD.send(None, dest=0, tag=WORK_REQ_TAG)
        ligandName = MPI.COMM_WORLD.recv(source=0, tag=MPI.ANY_TAG, status=wStatus)

    if wStatus.Get_tag() == TERMINATE_TAG:
        print "Worker {0} has terminated.\n".format(workerID)
    else:
        print "Worker {0} has received invalid Tag\n".format(workerID)

if __name__ == '__main__':
    main()