
function disp_pt(Ne)
    fps = 10
    Tf = 75

    avg_disp_mat = zeros(35,35)
    best_phase_mat = zeros(35,35)
    phases = [pt for pt = -pi:0.05:pi]
    disp_pt = zeros(size(phases))

    #set amount simulated
    global num_parts;
    num_parts = 3*Ne
    num_per = 3
    global dt;
    dt = 2/fps

    #find neighbors for every vertex
    global every_neighbor;
    every_neighbor = []
    for i=1:Nv
        # pos neighboring edges 
        push!(every_neighbor, [in_neighbors[i] ; out_neighbors[i]])
    end

    #set initial conditions
    global pos_0;
    pos_0 = []
    for i = 1:num_per
        append!(pos_0,i.*Len/(num_per+1))
    end

    global ind_0;
    ind_0 = []
    for i = 1:num_per
        append!(ind_0,[i for i = 1:Ne])
    end

    #num_samples = size([i for i = 0:1/fps:Tf],1)
    # Qin_t = [Vector{Float64}(undef,Ne) for _ in 1:num_samples]
    # Qout_t = [Vector{Float64}(undef,Ne) for _ in 1:num_samples]
    # A_t = [Vector{Float64}(undef,Ne) for _ in 1:num_samples]

    num_steps = 1+length(2:2:(length(A_t)-2))
    pos_time = zeros(Float64, num_steps, num_parts)
    ind_time = zeros(Int, num_steps, num_parts)

    pos_time[1,:] .= pos_0
    ind_time[1,:] .= ind_0

    Qmid = zeros(Ne,1); 

    mode_one = [0 for i = 0:1/fps:Tf]
    mode_two = [0 for i = 0:1/fps:Tf]
    radii_time_modes = zeros(size(mode_one,1),Ne)
    dV_modes = zeros(size(radii_time_modes,1),Ne)

    for x = 1:20
        println(x)
        for y = 1:20
            if x >= y
                #find maximum
                for z = 1:size(phases,1)
                    modes_init_conds(Qin_t, Qout_t, A_t, Ne-x, Ne-y, phases[z], fps, Tf, Qmid, mode_one, mode_two, radii_time_modes, dV_modes)
                    simulate(pos_time, ind_time, A_t, Qin_t, Qout_t, num_parts)
                    disp_pt[z] = find_average_displacement(pos_time, ind_time)
                end
                best_phase_mat[x,y] = phases[argmax(disp_pt)]
                avg_disp_mat[x,y] = maximum(disp_pt)
            end
        end
    end

    output = string("_DISPPT.txt")
    open(output,"w") do io
        writedlm(io,avg_disp_mat)
    end

    output = string("_PHASEPT.txt")
    open(output,"w") do io
        writedlm(io,best_phase_mat)
    end
    return disp
end

file_suffix = "sample_spring_one"
#p = [0.38 0.8]

global g; global avg_radii; global Len
g, avg_radii, Len = import_MATLAB_network(file_suffix, true)
params = graph_precomputation(g, avg_radii, Len)

global Ne;  global Nv; global avg_radii; global Len; global V0; global kappa; global α; global D; global β; global Deff; global Dt; global τ; global C0; global E; global B; global K; global k; global L; global Ldag; global Mq_inv; global Minv1; global sum_Minv1; global elastic_modes; global elastic_eigvals; global edgelist; global out_neighbors; global in_neighbors; global Ein; global Eout; global EinT; global EoutT; global surface_area; global radius_sq; global Char_L; global Char_t; global R0; global Cq; global mu; global b 
Ne, Nv, avg_radii, Len, V0, kappa, α, D, β, Deff, Dt, τ, C0, E, B, K, k, L, Ldag, Mq_inv, Minv1, sum_Minv1, elastic_modes, elastic_eigvals, edgelist, out_neighbors, in_neighbors, Ein, Eout, EinT, EoutT, surface_area, radius_sq, Char_L, Char_t, R0, Cq, mu, b = params


@time disp_max = disp_pt(Ne)
