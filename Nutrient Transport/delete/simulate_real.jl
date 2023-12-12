#Simulate nutrient transport
function real_init_conds(file_suffix, fps, sampleTime)

    global g; global avg_radii; global Len
    g, avg_radii, Len = import_MATLAB_network(file_suffix, true)
    params = graph_precomputation(g, avg_radii, Len)

    global Ne;  global Nv; global avg_radii; global Len; global V0; global kappa; global α; global D; global β; global Deff; global Dt; global τ; global C0; global E; global B; global K; global k; global L; global Ldag; global Mq_inv; global Minv1; global sum_Minv1; global elastic_modes; global elastic_eigvals; global edgelist; global out_neighbors; global in_neighbors; global Ein; global Eout; global EinT; global EoutT; global surface_area; global radius_sq; global Char_L; global Char_t; global R0; global Cq; global mu; global b 
    Ne, Nv, avg_radii, Len, V0, kappa, α, D, β, Deff, Dt, τ, C0, E, B, K, k, L, Ldag, Mq_inv, Minv1, sum_Minv1, elastic_modes, elastic_eigvals, edgelist, out_neighbors, in_neighbors, Ein, Eout, EinT, EoutT, surface_area, radius_sq, Char_L, Char_t, R0, Cq, mu, b = params


     #set amount simulated
     num_parts = 3*Ne
     num_per = 3
     dt = 2/fps

     #set initial conditions
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

    #get needed values
    # real 
    source = string("Raw data/",file_suffix,"/",file_suffix)
    destination = string("Processed Data/",file_suffix,"/",file_suffix)

    filename_radii = string(source,"_radii_time.txt")
    radii_data = readdlm(filename_radii) #.* 2 * 1.83 * 10^(-6)
    radii_time_real = zeros(length(split(radii_data[3],",")),Ne)

    for i = 1:Ne
        broken_up = split(radii_data[i],",")
        radii_data_curr = parse.(Float64,broken_up)
        radii_time_real[:,i] = radii_data_curr * 1.83 * 10^(-6) * (1/Char_L)
    end

    dV_real = zeros(size(radii_time_real,1),Ne)
    #real
    for i = 1:Ne
        #find flows 
        radii_time_curr = radii_time_real[:,i]

        #find dV
        volume_time = pi .* radii_time_curr.^2 .* Len[i]
        dV_real[2:end-1,i] = (volume_time[3:end] - volume_time[1:end-2])/(2*sampleTime/60)
        dV_real[1,i] = (volume_time[2] - volume_time[1]) / (sampleTime/60)
        dV_real[end,i] = (volume_time[end] - volume_time[end-1]) / (sampleTime/60)
    end

    Qin_t = [] 
    Qout_t = []
    A_t = []
    for i = 1:size(radii_time_real,1)
        Qmid =  mid_flows(dV_real[i,:])
        Qin, Qout = in_out_flows(Qmid, dV_real[i,:])
        push!(Qin_t,Qin)
        push!(Qout_t,Qout)
        push!(A_t,pi .* radii_time_real[i,:].^2)
    end
    return num_parts, num_per, dt, pos_0, ind_0, every_neighbor, Qin_t, Qout_t, A_t
end