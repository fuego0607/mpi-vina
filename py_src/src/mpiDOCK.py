''' File      : mpiDOCK.py
*   Author    : Jeremy Hofer and Abhay Aradhya
*   Purpose   : Python MPI based parallel master/slave manager of Autodock Vina and various Vina rescoring functions
'''
from mpi4py import MPI
from collections import deque
from glob import glob
from os import system, access, X_OK
from datetime import datetime

MASTER = 0

COMPUTE_TAG = 11
TERMINATE_TAG = 22
WORK_REQ_TAG = 33

config_file = "mpidock.config"

def abort_mpi(error_message):
    print error_message
    print "\nTerminating mpiDOCK with error code 1."
    MPI.COMM_WORLD.Abort(1)
    exit()

def main():
    numProcs = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()

    if numProcs < 2:
        abort_mpi("Not enough processors! You must use at least 2 processors.")

    #Only master processor will read the ligandlist file and will make the work pool.
    if rank == MASTER:
        print "Master processor initializing mpiDOCK.".format(MASTER)
        print "\nReading configuration file..."

        configuration = {"vina": None, "vina_config": None, "run_vina": False, "receptors_dir": None, "ligands_dir": None,
                         "job_name": datetime.now().strftime("%Y-%m-%d_%H%M%S"), "ligand_db_name": None, "output_dir": None,
                         "clean_temp_files": None}
        
        #try to open and read in configuration file. abort MPI if errors occur
        try:
            file = open(config_file, "r")
        except (FileNotFoundError, IOError) as e:
            abort_mpi("Cannot find or open the configuration file. Please make sure mpidock.config is in the same directory as mpiDOCK.py.")
        else:
            lines = [x.strip().split("=") for x in file.readlines() if x[0] is not "#"]
            if len(lines) is 0:
                abort_mpi("Error reading configuration file. No configuration parameters seen. Make sure they are not commented out with '#'.")

            for line in lines:
                if line[0] in configuration:
                    if line[0] in ["run_vina", "clean_temp_files"]:
                        if line[1] == "yes":
                            line[1] = True
                        else:
                            line[1] = False

                    configuration[line[0]] = line[1]
                else:
                    abort_mpi("Invalid key/variable {0} in configuration file. Please review valid keys/variables and try again .")

        #print configuration information
        print "Configuration read OK. Configured with the following parameters:"

        for key, value in configuration.iteritems():
            print "{0} : {1}".format(key, value)

        #perform checks to ensure all necessary files are present or will be present given the set of configured actions to perform (dock, rescore, etc.)
        print "\nVerifying all required files are present..."

        if configuration['run_vina']:
            #check for vina executable
            exe = access(configuration['vina'], X_OK)

            if not exe:
                abort_mpi("Issue validating Vina executable. File does not exist at location or does not have executable permissions.")

            #check for configuration file
            try:
                open(configuration['vina_config'], "r")
            except (FileNotFoundError, IOError) as e:
                abort_mpi("Vina configuration file not found at given path or does not have read permissions.")

            print "Vina files verified."

        print "\nGenerating receptor and ligand pools..."
        receptors = glob(configuration['receptors_dir']+"/*.pdbqt")
        ligands = glob(configuration['ligands_dir']+"/*.pdbqt")
        num_receptors = len(receptors)
        num_ligands = len(ligands)

        if len(receptors) is 0 or len(ligands) is 0:
            abort_mpi("Cannot operate. {0} receptors found and {1} ligands found. Must have at least 1 receptor and 1 ligand.".format(num_receptors, num_ligands))
        else:
            print "{0} receptors found and {1} ligands found.".format(num_receptors, num_ligands)

        #initialize job queue
        queue = deque([])

        #append items to the queue
        for i in receptors:
            for j in ligands:
                queue.append({"receptor": i, "ligand": j})

        total_jobs = len(queue)

        print "\n---- STARTING mpiDOCK job {0} ----".format(configuration['job_name'])

    MPI.COMM_WORLD.Barrier()

    if rank == MASTER:
        startTime = MPI.Wtime(); #start timer.
        mpiVinaManager(numProcs, queue, configuration);   #Master processor will play the role of mpiVINA manager.
    else:
        mpiVinaWorker(rank, configuration);    #All other processors will play the role of mpiVINA worker.

    MPI.COMM_WORLD.Barrier()

    if rank == MASTER:
        endTime = MPI.Wtime(); #end timer.
        print "\n.........................................." 
        print "   Number of workers       = {0}".format(numProcs - 1)
        print "   Number of dockings        = {0}".format(total_jobs)
        print "   Total time required     = {0} seconds.".format(endTime - startTime)
        print ".........................................."

    MPI.Finalize()

def mpiVinaManager(numProcs, queue, configuration):
    print "Config in manager is {1}".format(configuration)
    mStatus = MPI.Status()
    while len(queue) > 0:
        MPI.COMM_WORLD.recv(source=MPI.ANY_SOURCE, tag=WORK_REQ_TAG, status=mStatus)
        print "Worker {0} requesting work\n".format(mStatus.Get_source())
        MPI.COMM_WORLD.send(queue.popleft(), dest=mStatus.Get_source(), tag=COMPUTE_TAG)

    for i in range(numProcs - 1):
        MPI.COMM_WORLD.recv(source=MPI.ANY_SOURCE, tag=WORK_REQ_TAG, status=mStatus)
        print "Sending termination to worker {0}\n".format(mStatus.Get_source())
        MPI.COMM_WORLD.send(None, dest=mStatus.Get_source(), tag=TERMINATE_TAG)

def mpiVinaWorker(workerID, configuration):
    print "Worker {0} has started.\n".format(workerID)
    wStatus = MPI.Status()

    MPI.COMM_WORLD.send(None, dest=0, tag=WORK_REQ_TAG)
    job_info = MPI.COMM_WORLD.recv(source=0, tag=MPI.ANY_TAG, status=wStatus)

    while wStatus.Get_tag() == COMPUTE_TAG:
        print "Worker = {0} : receptor {1} ligand {2} is processing...\n".format(workerID, job_info['receptor'], job_info['ligand'])
        #insert vina command
        MPI.COMM_WORLD.send(None, dest=0, tag=WORK_REQ_TAG)
        ligandName = MPI.COMM_WORLD.recv(source=0, tag=MPI.ANY_TAG, status=wStatus)

    if wStatus.Get_tag() == TERMINATE_TAG:
        print "Worker {0} has terminated.\n".format(workerID)
    else:
        print "Worker {0} has received invalid Tag\n".format(workerID)

if __name__ == '__main__':
    main()