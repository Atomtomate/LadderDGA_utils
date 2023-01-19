using HDF5

file = ARGS[1]

f = h5open(file, "r")

hist = read(f["dmft-last/ineq-001/hist/value"])
val, ind = findmax(hist)
max_N = ind[1]
max_accept_rem = maximum(read(f["dmft-last/ineq-001/accept-rem/value"]))
N_corr = ceil(Int,max_N/max_accept_rem)
println(N_corr)
