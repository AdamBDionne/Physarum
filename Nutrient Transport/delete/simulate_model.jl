#Simulate nutrient transport
include("fxns/helper_functions.jl")
include("fxns/model_helper_functions.jl")


function model_init_conds(Ne, fps)
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
    #model 
    Qin_t = []; Qout_t = [];
    for x = 1:length(sol.u)
        Qin_tc, Qout_tc = findFlows(sol.u[x])
        push!(Qin_t, Qin_tc); push!(Qout_t, Qout_tc)
    end
   # Qin_t, Qout_t = findFlows.(sol.u)
    #Qin_t = findFlowsOut.(sol.u)
    #Qout_t = findFlowsIn.(sol.u)
    A_t = []
    for i = 1:length(sol.t)
        push!(A_t,sol.u[i][1:Ne] ./ (Len))
    end
    return A_t, Qin_t, Qout_t, pos_0, ind_0, num_parts, num_per, every_neighbor, dt, pos_0, ind_0
end

#@time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, pos_0, ind_0)
