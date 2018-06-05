''' File      : mpiDOCK.py
*   Author    : Jeremy Hofer and Abhay Aradhya
*   Purpose   : Python MPI based parallel master/slave manager of Autodock Vina and various Vina rescoring functions
'''
from mpi4py import MPI
from collections import deque
from glob import glob
from os import system, access, X_OK, R_OK, makedirs
from os.path import exists, basename, splitext
from datetime import datetime, timedelta
from subprocess import call
import score_commands as commands
import sys, getopt

MASTER = 0

COMPUTE_TAG = 11
TERMINATE_TAG = 22
WORK_REQ_TAG = 33
DOCK_FAIL_TAG = 44

config_file = "config/mpidock.config"

def abort_mpi(error_message):
    print error_message
    print "\nTerminating mpiDOCK with error code 1."
    MPI.COMM_WORLD.Abort(1)
    exit()

def main(argv):
    numProcs = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    configuration = None

    if numProcs < 2:
        abort_mpi("Not enough processors! You must use at least 2 processors.")

    #Only master processor will read the ligandlist file and will make the work pool.
    if rank == MASTER:
        print "Master processor initializing mpiDOCK.".format(MASTER)
        print "\nReading configuration file..."

        configuration = {"vina": None, "vina_config": None, "run_vina": False, "receptors_dir": None, "ligands_dir": None,
                         "job_name": datetime.now().strftime("%Y-%m-%d_%H%M%S_%f"), "ligand_db_name": None, "output_dir": None,
                         "clean_temp_files": None, "run_dsx": False, "run_xscore": False, "run_nnscore": False, "run_rfscore": False,
                         "pdb_potentials_dir": None, "xscore_base_config": None, "networks_dir": None, "dsx": None, "xscore": None,
                         "nnscore": None, "rfscore": None, "job_out_dir": None, "temp_out_dir": None, "vina_out_dir": None,
                         "dsx_out_dir": None, "xscore_out_dir": None, "fixmol2_dir": None, "fixpdb_dir": None, "mol2_out_dir": None,
                         "nnscore_dir": None, "rfscore_dir": None, "obabel": None}
        
        #try to open and read in configuration file. abort MPI if errors occur
        try:
            file = open(config_file, "r")
        except IOError as e:
            abort_mpi("Cannot find or open the configuration file. Please make sure mpidock.config is in the same directory as mpiDOCK.py.")
        else:
            lines = [x.strip().split("=") for x in file.readlines() if x[0] is not "#"]
            if len(lines) is 0:
                abort_mpi("Error reading configuration file. No configuration parameters seen. Make sure they are not commented out with '#'.")

            for line in lines:
                if line[0] in configuration:
                    if line[0] in ["run_vina", "clean_temp_files", "run_dsx", "run_xscore", "run_nnscore", "run_rfscore"]:
                        if line[1] == "yes":
                            line[1] = True
                        else:
                            line[1] = False

                    configuration[line[0]] = line[1]
                else:
                    abort_mpi("Invalid key/variable {0} in configuration file. Please review valid keys/variables and try again.".format(line[0]))

        #print configuration information
        print "Configuration read OK. Configured with the following parameters:"

        for key, value in configuration.iteritems():
            print "{0} : {1}".format(key, value)

        #perform checks to ensure all necessary files are present or will be present given the set of configured actions to perform (dock, rescore, etc.)
        print "\nVerifying all required files and directories are present..."

        configuration['job_out_dir'] = configuration['output_dir']+"/"+configuration['job_name']
        configuration['temp_out_dir'] = "./temp/"+configuration['job_name']

        if configuration['run_vina']:
            #check for vina executable
            exe = access(configuration['vina'], X_OK)

            if not exe:
                abort_mpi("Issue validating Vina executable. File does not exist at location or does not have executable permissions.")

            #check for configuration file
            try:
                f = open(configuration['vina_config'], "r")
                f.close()
            except IOError as e:
                abort_mpi("Vina configuration file not found at given path or does not have read permissions.")

            configuration['vina_out_dir'] = configuration['job_out_dir']+"/vina_out"

            #make vina output directory
            if not exists(configuration['vina_out_dir']):
                makedirs(configuration['vina_out_dir'])

            print "Vina files verified."

        if configuration['run_dsx']:
            #check for dsx executable
            exe = access(configuration['dsx'], X_OK)

            if not exe:
                abort_mpi("Issue validating DSX executable. File does not exist at location or does not have executable permissions.")

            #check for populated potentials directory
            if len(glob(configuration['pdb_potentials_dir']+"/*")) is 0:
                abort_mpi("Issue with pdb potentials directory. Path is not valid or there are no files present.")

            configuration['dsx_out_dir'] = configuration['job_out_dir']+"/dsx_out"

            if not exists(configuration['dsx_out_dir']):
                makedirs(configuration['dsx_out_dir'])

            print "DSX files verified."

        if configuration['run_xscore']:
            #check for vina executable
            exe = access(configuration['xscore'], X_OK)

            if not exe:
                abort_mpi("Issue validating XScore executable. File does not exist at location or does not have executable permissions.")

            configuration['xscore_out_dir'] = configuration['job_out_dir']+"/xscore_out"
            configuration['fixmol2_dir'] = configuration['temp_out_dir']+"/fixmol2"
            configuration['fixpdb_dir'] = configuration['temp_out_dir']+"/fixpdb"

            if not exists(configuration['xscore_out_dir']):
                makedirs(configuration['xscore_out_dir'])

            if not exists(configuration['fixmol2_dir']):
                makedirs(configuration['fixmol2_dir'])

            if not exists(configuration['fixpdb_dir']):
                makedirs(configuration['fixpdb_dir'])

            print "XScore files verified."

        if configuration['run_dsx'] or configuration['run_xscore']:
            configuration['mol2_out_dir'] = configuration['temp_out_dir']+"/mol2_out"

            if not exists(configuration['mol2_out_dir']):
                makedirs(configuration['mol2_out_dir'])

        if configuration['run_nnscore']:
            #check for dsx executable
            exe = access(configuration['nnscore'], R_OK)

            if not exe:
                abort_mpi("Issue validating NNScore script. File does not exist at location or does not have read permissions.")

            #check for populated potentials directory
            if len(glob(configuration['networks_dir']+"/*")) is 0:
                abort_mpi("Issue with neural networks directory. Path is not valid or there are no files present.")

            configuration['nnscore_dir'] = configuration['job_out_dir']+"/nnscore_out"

            if not exists(configuration['nnscore_dir']):
                makedirs(configuration['nnscore_dir'])

            print "NNScore files verified."

        if configuration['run_rfscore']:
            #check for vina executable
            exe = access(configuration['rfscore'], X_OK)

            if not exe:
                abort_mpi("Issue validating RFScore executable. File does not exist at location or does not have executable permissions.")

            configuration['rfscore_dir'] = configuration['job_out_dir']+"/rfscore_out"

            if not exists(configuration['rfscore_dir']):
                makedirs(configuration['rfscore_dir'])

            print "RFScore files verified."

        print "\nGenerating receptor and ligand pools..."
        receptors = glob(configuration['receptors_dir']+"/*.pdbqt")
        receptors_pdb = glob(configuration['receptors_dir']+"/*.pdb")
        ligands = glob(configuration['ligands_dir']+"/*.pdbqt")
        num_receptors = len(receptors)
        num_ligands = len(ligands)

        if len(receptors) is 0 or len(ligands) is 0:
            abort_mpi("Cannot operate vina. {0} receptors found and {1} ligands found. Must have at least 1 receptor and 1 ligand.".format(num_receptors, num_ligands))
        elif len(receptors_pdb) is 0 and (configuration['run_dsx'] or configuration['run_xscore']):
            abort_mpi("Cannot operate dsx or xscore. No receptor files found in pdb format")
        else:
            print "{0} receptors found and {1} ligands found.".format(num_receptors, num_ligands)

        #initialize job queue
        queue = deque([])

        #append items to the queue
        for i in receptors:
            if configuration['run_dsx'] or configuration['run_xscore']:
                #ensure matching pdb is found
                receptor_pdb = i.replace(".pdbqt", ".pdb")
                if receptor_pdb not in receptors_pdb:
                    abort_mpi("No matching .PDB receptor file found for {0}").format(i)
            for j in ligands:
                queue.append({"receptor": i, "ligand": j, "receptor_pdb": receptor_pdb})

        total_jobs = len(queue)
    
        print "\n---- STARTING mpiDOCK job {0} ----".format(configuration['job_name'])

    configuration = MPI.COMM_WORLD.bcast(configuration, 0)
    MPI.COMM_WORLD.Barrier()

    if rank == MASTER:
        startTime = MPI.Wtime(); #start timer.
        #mpiVinaManager(numProcs - 1, queue, configuration);   #Master processor will play the role of mpiVINA manager.
        pass
    else:
        #mpiVinaWorker(rank, configuration);    #All other processors will play the role of mpiVINA worker.
        pass

    MPI.COMM_WORLD.Barrier()

    if rank == MASTER:
        endTime = MPI.Wtime(); #end timer.
        print "\n.........................................." 
        print "   Number of workers       = {0}".format(numProcs - 1)
        print "   Number of dockings        = {0}".format(total_jobs)
        print "   Total time required     = {0}".format(str(timedelta(seconds=(endTime - startTime))))
        print ".........................................."

    MPI.Finalize()

def mpiVinaManager(numWorkers, queue, configuration):
    mStatus = MPI.Status()
    fail_list = []
    while len(queue) > 0 and numWorkers > 0:
        job_done = MPI.COMM_WORLD.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=mStatus)

        if mStatus.Get_tag() == DOCK_FAIL_TAG:
            fail_list.append(job_done)
            MPI.COMM_WORLD.send(queue.popleft(), dest=mStatus.Get_source(), tag=COMPUTE_TAG)
        elif mStatus.Get_tag() == WORK_REQ_TAG:
            #print "Worker {0} requesting work\n".format(mStatus.Get_source())
            MPI.COMM_WORLD.send(queue.popleft(), dest=mStatus.Get_source(), tag=COMPUTE_TAG)
        else:
            #print "Unknown tag {0} from {1}. Terminating process.".format(mStatus.Get_tag(), mStatus.Get_source())
            MPI.COMM_WORLD.send(None, dest=mStatus.Get_source(), tag=TERMINATE_TAG)
            numWorkers = numWorkers - 1
    
    while numWorkers > 0:
        job_done = MPI.COMM_WORLD.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=mStatus)

        if mStatus.Get_tag() == DOCK_FAIL_TAG:
            fail_list.append(job_done)

        #print "Sending termination to worker {0}\n".format(mStatus.Get_source())
        MPI.COMM_WORLD.send(None, dest=mStatus.Get_source(), tag=TERMINATE_TAG)
        numWorkers = numWorkers - 1

def mpiVinaWorker(workerID, configuration):
    #print "Worker {0} has started.\n".format(workerID)
    #source the user's bash_profile on each node to ensure all environment variables are properly set
    wStatus = MPI.Status()

    MPI.COMM_WORLD.send(None, dest=0, tag=WORK_REQ_TAG)
    job_info = MPI.COMM_WORLD.recv(source=0, tag=MPI.ANY_TAG, status=wStatus)

    while wStatus.Get_tag() == COMPUTE_TAG:
        #print "Worker = {0} : receptor {1} ligand {2} is processing...\n".format(workerID, job_info['receptor'], job_info['ligand'])

        receptor = basename(splitext(job_info['receptor'])[0])
        ligand = basename(splitext(job_info['ligand'])[0])
        output_base = receptor + "_" + ligand

        vina_out = configuration['vina_out_dir']+"/"+output_base+".pdbqt"

        mol2_out = configuration['mol2_out_dir']+"/"+output_base+".mol2"

        dsx_out = configuration['dsx_out_dir']+"/"+output_base+".dsx_out"
        
        xscore_out = configuration['xscore_out_dir']+"/"+output_base+".xscore_out"
        
        fixmol2_out = configuration['fixmol2_dir']+"/"+output_base+".mol2"
        
        fixpdb_out = configuration['fixpdb_dir']+"/"+output_base+".pdb"

        nnscore_out = configuration['nnscore_dir']+"/"+output_base+".nnscore_out"

        rfscore_out = configuration['rfscore_dir']+"/"+output_base+".rfscore_out"

        if configuration['run_vina']:
            commands.runVina(configuration['vina'], configuration['vina_config'], job_info['receptor'], job_info['ligand'], vina_out)

        if configuration['run_dsx'] or configuration['run_xscore']:
            #convert vina_out to mol2_out
            if exists(vina_out):
                commands.convertToMol2(configuration['obabel'], vina_out, mol2_out)

                if configuration['run_dsx']:
                    
                    commands.runDSXScore(configuration['dsx'], job_info['receptor_pdb'], mol2_out, configuration['pdb_potentials_dir'], dsx_out)

                if configuration['run_xscore']:
                    
                    commands.runXScore(configuration['xscore'], job_info['receptor_pdb'], fixpdb_out, mol2_out, fixmol2_out, xscore_out)

        if configuration['run_nnscore']:
            if exists(vina_out):
                commands.runNNScore(configuration['nnscore'], job_info['receptor'], vina_out, configuration['networks_dir'], nnscore_out)


        if configuration['run_rfscore']:
            if exists(vina_out):
                commands.runRFScore(configuration['rfscore'], job_info['receptor'], vina_out, rfscore_out)

        MPI.COMM_WORLD.send(None, dest=0, tag=WORK_REQ_TAG)
        job_info = MPI.COMM_WORLD.recv(source=0, tag=MPI.ANY_TAG, status=wStatus)

    if wStatus.Get_tag() == TERMINATE_TAG:
        #print "Worker {0} has terminated.\n".format(workerID)
        pass
    else:
        #print "Worker {0} has received invalid Tag\n".format(workerID)
        pass

if __name__ == '__main__':
    main(sys.argv[1:])
