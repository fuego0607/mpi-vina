from mpi4py import MPI
a = [1,2,3]
if MPI.COMM_WORLD.rank == 0:
        MPI.COMM_WORLD.send(a, dest = 1)
else:
        a = MPI.COMM_WORLD.recv(source = 0)
        print a

        
from mpi4py import MPI
import sys

size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()

sys.stdout.write(
    "Hello, World! I am process %d of %d on %s.\n"
    % (rank, size, name))
