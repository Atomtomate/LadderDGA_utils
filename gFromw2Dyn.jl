using HDF5

f = h5open(ARGS[1])

g_in = read(f["dmft-last/ineq-001/giw/value"])
