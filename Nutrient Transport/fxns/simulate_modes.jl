#setting initial conditions to simulate nutrient transport for two modes with a given phase shift phi
function modes_init_conds(mode_num_1, mode_num_2, phi, fps, Tf)
    #init 
    num_samples = size([i for i = 0:1/fps:Tf],1)
    Qin_t = [Vector{Float64}(undef,Ne) for _ in 1:num_samples]
    Qout_t = [Vector{Float64}(undef,Ne) for _ in 1:num_samples]
    A_t = [Vector{Float64}(undef,Ne) for _ in 1:num_samples]
    dV_modes = zeros(size(radii_time_modes,1),Ne)


    mode_one =  [(0.12 .* elastic_modes[:,mode_num_1] .* sin((2*pi*i/8) + phi)) for i = 0:1/fps:Tf]
    mode_two =  [(0.12 .* elastic_modes[:,mode_num_2] .* sin(2*pi*i/8)) for i = 0:1/fps:Tf]
    
    radii_time_modes = zeros(size(mode_one,1),Ne)
    for i = 1:size(mode_one,1)
        radii_time_modes[i,:] = avg_radii .* (1 .+ mode_one[i] .+ mode_two[i])
    end

    for i = 1:Ne
        #find flows 
        radii_time_curr = radii_time_modes[:,i]

        #find dV
        volume_time = pi .* radii_time_curr.^2 .* Len[i]
        dV_modes[2:end-1,i] = (volume_time[3:end] - volume_time[1:end-2])/(2/fps)
        dV_modes[1,i] = (volume_time[2] - volume_time[1]) / (1/fps)
        dV_modes[end,i] = (volume_time[end] - volume_time[end-1]) / (1/fps)
    end

    
    for i = 1:size(radii_time_modes,1)
        Qmid = mid_flows(dV_modes[i,:])
        Qin_t[i] = Qmid .+ (0.5 .* dV_modes[i,:])
        Qout_t[i] = Qmid .- (0.5 .* dV_modes[i,:])
        A_t[i] = pi .* radii_time_modes[i,:].^2
    end



    num_per = 3
    pos_0 = []

    for i = 1:num_per
        append!(pos_0,i.*Len/(num_per+1))
    end
    ind_0 = []
    for i = 1:num_per
        append!(ind_0,[i for i = 1:Ne])
    end
    #find neighbors for every vertex
    every_neighbor = []
    for i=1:Nv
        # pos neighboring edges 
        push!(every_neighbor, [in_neighbors[i] ; out_neighbors[i]])
    end

    return Qin_t, Qout_t, A_t, radii_time_modes, pos_0, ind_0, every_neighbor
end