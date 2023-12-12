using SimpleWeightedGraphs
using LinearAlgebra
using Statistics
using LightGraphs
using StatsBase



#find particles current velocity
function find_vel(A_t, Qin_t, Qout_t, t_ind, part_pos, part_ind, num_parts)
    A = A_t[t_ind]
    Qin = Qin_t[t_ind]
    Qout = Qout_t[t_ind]

    vel = zeros(num_parts)
    counter = 1
    for ind in part_ind
        Q = (Qout[ind] - Qin[ind])*(part_pos[counter]/Len[ind]) + Qin[ind]
        vel[counter] = 2*Q/A[ind]
        counter = counter+1
    end

    return vel
end

#update position of particle using RK4
function update_pos(A_t, Qin_t, Qout_t, t_ind, part_pos, part_ind, num_parts)
    k1 = find_vel(A_t, Qin_t, Qout_t, t_ind,part_pos, part_ind, num_parts) 
    k2 = find_vel(A_t, Qin_t, Qout_t, t_ind-1, part_pos+(dt .* k1 ./2), part_ind, num_parts) 
    k3 = find_vel(A_t, Qin_t, Qout_t, t_ind-1, part_pos+(dt .* k2 ./ 2), part_ind, num_parts)
    k4 = find_vel(A_t, Qin_t, Qout_t, t_ind+2, part_pos+(dt .*k3), part_ind, num_parts)
    return part_pos .+ (1/6)*dt*(k1+2k2+2k3+k4)
end

#simulate nutrient transport
function simulate(A_t, Qin_t, Qout_t, num_parts)
    num_steps = 1+length(2:2:(length(A_t)-2))
    pos_time = zeros(Float64, num_steps, num_parts)
    ind_time = zeros(Int, num_steps, num_parts)

    pos_time[1,:] .= pos_0
    ind_time[1,:] .= ind_0

    for t_ind = 2:2:(length(A_t)-2)
        s_ind = Int((t_ind/2)+1)
        pos_time[s_ind,:] = update_pos(A_t, Qin_t, Qout_t, t_ind, pos_time[s_ind-1,:], ind_time[s_ind-1,:], num_parts)
        counter = 1
        ind_time[s_ind,:] = ind_time[s_ind-1,:]
        for pos in pos_time[s_ind,:]
            if pos < 0 #right node, ie node (j): Qin
                #current vessel index 
                ind = ind_time[s_ind-1,counter]
                node = edgelist[ind,2]

                #possible next indices 
                poss_inds = copy(every_neighbor[node])
                deleteat!(poss_inds, poss_inds .== ind);

                #find flow from each edge exiting our starting node 
                poss_flows = Float64[]
                valid_poss_inds = Int64[]
                for curr_ind in poss_inds
                    nodes = edgelist[curr_ind,:]
                    if nodes[1] == node
                        #must be flowing out of starting node 
                        if Qout_t[t_ind][curr_ind] <= 0 
                            push!(poss_flows, abs(Qout_t[t_ind][curr_ind]))
                            push!(valid_poss_inds, curr_ind)
                        end
                    else
                        if Qin_t[t_ind][curr_ind] >= 0 
                            push!(poss_flows, abs(Qin_t[t_ind][curr_ind]))
                            push!(valid_poss_inds, curr_ind)
                        end
                    end
                end
                if isempty(poss_flows)
                    #no valid place to flow 
                    #stay in same vessel, go a little farther up 
                    pos_time[s_ind,counter] = 1e-5
                    #integrator.u[idx] = 1e-5#Len[ind]/100
                else
                    #weighted probability 
                    new_ind = StatsBase.sample(valid_poss_inds, Weights(poss_flows))
                    ind_time[s_ind,counter] = new_ind
                    
                    # are we at the end or at the beginning?
                    if edgelist[new_ind,1] == node
                        pos_time[s_ind,counter] = Len[new_ind] - 1e-5
                    else
                        pos_time[s_ind,counter] = 0.0 + 1e-5
                    end
                end
            end
            if pos > Len[ind_time[s_ind-1,counter]]
                #left node, ie node (i): Qout
                ind = ind_time[s_ind-1,counter]
                node = edgelist[ind,1]

                #possible next indices 
                poss_inds = copy(every_neighbor[node])
                deleteat!(poss_inds, poss_inds .== ind);

                #find flow from each edge exiting our starting node 
                poss_flows = Float64[]
                valid_poss_inds = Int64[]
                for curr_ind in poss_inds
                    nodes = edgelist[curr_ind,:]
                    if nodes[1] == node
                        #must be flowing out of starting node 
                        if Qout_t[t_ind][curr_ind] <= 0 
                            push!(poss_flows, abs(Qout_t[t_ind][curr_ind]))
                            push!(valid_poss_inds, curr_ind)
                        end
                    else
                        if Qin_t[t_ind][curr_ind] >= 0 
                            push!(poss_flows, abs(Qin_t[t_ind][curr_ind]))
                            push!(valid_poss_inds, curr_ind)
                        end
                    end
                end
                if isempty(poss_flows)
                    #no valid place to flow 

                    #stay in same vessel, go a little farther down
                    pos_time[s_ind,counter] = Len[ind] - 1e-5#(Len[ind]/100)
                else
                    #weighted probability 
                    new_ind = StatsBase.sample(valid_poss_inds, Weights(poss_flows))
                    ind_time[s_ind,counter] = new_ind

                    # are we at the end or at the beginning?
                    if edgelist[new_ind,1] == node
                        pos_time[s_ind,counter] = Len[new_ind] - 1e-5#1e-8
                    else
                        pos_time[s_ind,counter] = 0.0 + 1e-5#1e-8
                    end
                end
            end
            counter = counter + 1
        end
    end
    return pos_time, ind_time
end


#find flows given current volumes and concentrations
function findFlows(u)
    V = u[1:Ne] #volume 
    C = u[Ne+1:2Ne] #solute amount
    radius_sq = (V) ./ (pi .* Len)
    R = sqrt.(radius_sq)
    ϵ = (R .- R0)./R0 #percent change from equilibrium volume

    #Volume dynamics
                  #  restoring force      nonlinear fudge    actomyosin contractions 
    r = Mq_inv \ (-(k .*ϵ) .- (kappa .*(ϵ).^3) .- ( α .*(C./C0).*(1 .- ϵ/ϵ_s)))
    μ = sum(r) / sum_Minv1
    Vdot = r .- (Minv1 .* μ)

    #Flows 
    Qmid = mid_flows(Vdot)
    return in_out_flows(Qmid, Vdot)
end


# #find flows given current volumes and concentrations
function findFlowsOut(u)
    V = u[1:Ne] #volume 
    C = u[Ne+1:2Ne] #solute amount
    radius_sq = (V) ./ (pi .* Len)
    R = sqrt.(radius_sq)
    ϵ = (R .- R0)./R0 #percent change from equilibrium volume

    #Volume dynamics
                  #  restoring force      nonlinear fudge    actomyosin contractions 
    r = Mq_inv \ (-(k .*ϵ) .- (kappa .*(ϵ).^3) .- ( α .*(C./C0).*(1 .- ϵ/ϵ_s)))
    μ = sum(r) / sum_Minv1
    Vdot = r .- (Minv1 .* μ)

    #Flows 
    Qmid = mid_flows(Vdot)
    _, Qout = in_out_flows(Qmid, Vdot)
    return Qout
end

