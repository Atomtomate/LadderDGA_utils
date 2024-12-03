using HDF5


β, g4iw, giw, vert_upup, vert_updo, ω_arr, ν_arr = h5open("Vertex_Z_new.hdf5", "r") do f
    β       = read_attribute(f[".config"], "general.beta")
    g4iw    = read(f["stat-last/ineq-001/g4iw/value"])
    giw_in  = read(f["stat-last/ineq-001/giw/value"])
    giw  = (giw_in[:,1,1] .+ giw_in[:,2,1]) ./ 2
    ω_arr = read(f[".axes/iwb-g4"])
    ν_arr = read(f[".axes/iwf-g4"])
    n0 = floor(Int,length(giw)/2)+1
    vert_upup = β .* (g4iw[:,:,:,1,1,1,1] .+ g4iw[:,:,:,2,1,2,1]) ./ 2
    vert_updo = β .* (g4iw[:,:,:,1,1,2,1] .+ g4iw[:,:,:,2,1,1,1]) ./ 2
    β, g4iw, giw, vert_upup, vert_updo, ω_arr, ν_arr
end
