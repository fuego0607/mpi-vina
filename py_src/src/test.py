from mpi4py import MPI
a = [1,2,3]
if MPI.COMM_WORLD.rank == 0:
        MPI.COMM_WORLD.send(a, dest = 1)
else:
        a = MPI.COMM_WORLD.recv(source = 0)
        print a