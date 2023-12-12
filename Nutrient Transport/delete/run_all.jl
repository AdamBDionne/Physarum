
# function disp_pt(Ne)
#     fps = 10
#     Tf = 75

#     avg_disp_mat = zeros(35,35)
#     best_phase_mat = zeros(35,35)
#     phases = [pt for pt = -pi:0.05:pi]
#     disp_pt = zeros(size(phases))

    #set amount simulated
    # global num_parts;
    # num_parts = 3*Ne
    # num_per = 3
    # global dt;
    # dt = 2/fps

    # #find neighbors for every vertex
    # global every_neighbor;
    # every_neighbor = []
    # for i=1:Nv
    #     # pos neighboring edges 
    #     push!(every_neighbor, [in_neighbors[i] ; out_neighbors[i]])
    # end


    # #set initial conditions
    # global pos_0;
    # pos_0 = []
    # for i = 1:num_per
    #     append!(pos_0,i.*Len/(num_per+1))
    # end

    # global ind_0;
    # ind_0 = []
    # for i = 1:num_per
    #     append!(ind_0,[i for i = 1:Ne])
    # end

    # #num_samples = size([i for i = 0:1/fps:Tf],1)
    # # Qin_t = [Vector{Float64}(undef,Ne) for _ in 1:num_samples]
    # # Qout_t = [Vector{Float64}(undef,Ne) for _ in 1:num_samples]
    # # A_t = [Vector{Float64}(undef,Ne) for _ in 1:num_samples]

    # num_steps = 1+length(2:2:(length(A_t)-2))
    # pos_time = zeros(Float64, num_steps, num_parts)
    # ind_time = zeros(Int, num_steps, num_parts)

    # pos_time[1,:] .= pos_0
    # ind_time[1,:] .= ind_0

#     Qmid = zeros(Ne,1); 

#     mode_one = [0 for i = 0:1/fps:Tf]
#     mode_two = [0 for i = 0:1/fps:Tf]
#     radii_time_modes = zeros(size(mode_one,1),Ne)
#     dV_modes = zeros(size(radii_time_modes,1),Ne)

#     x = 2; y = 1;
#     modes_init_conds(Qin_t, Qout_t, A_t, Ne-x, Ne-y, pi/2, fps, Tf, Qmid, mode_one, mode_two, radii_time_modes, dV_modes)
#     simulate(pos_time, ind_time, A_t, Qin_t, Qout_t, num_parts)
#     disp = find_average_displacement(pos_time, ind_time)

#     # for x = 1:1
#     #     println(x)
#     #     for y = 1:1
#     #         if x >= y
#     #             #find maximum
#     #             for z = 1:size(phases,1)
#     #                 x = 2 ; y = 1 ; 
#     #                 modes_init_conds(Qin_t, Qout_t, A_t, Ne-x, Ne-y, phases[z], fps, Tf, Qmid, mode_one, mode_two, radii_time_modes, dV_modes)
#     #                 simulate(pos_time, ind_time, A_t, Qin_t, Qout_t, num_parts)
#     #                 disp_pt[z] = find_average_displacement(pos_time, ind_time)
#     #             end
#     #             best_phase_mat[x,y] = phases[argmax(disp_pt)]
#     #             avg_disp_mat[x,y] = maximum(disp_pt)
#     #         end
#     #     end
#     # end

#     # output = string("_DISPPT.txt")
#     # open(output,"w") do io
#     #     writedlm(io,avg_disp_mat)
#     # end

#     # output = string("_PHASEPT.txt")
#     # open(output,"w") do io
#     #     writedlm(io,best_phase_mat)
#     # end
#     return disp
# end

# @time disp_max = disp_pt(Ne)


# #@time num_parts, num_per, dt, pos_0, ind_0, every_neighbor, Qin_t, Qout_t, A_t = real_init_conds(file_suffix, fps, sampleTime)
# num_parts, num_per, dt, pos_0, ind_0, every_neighbor, Qin_t, Qout_t, A_t = modes_init_conds(Ne-1, Ne-2, x_axis_model[argmax(t)], Ne, elastic_modes, fps, Tf)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, pos_0, ind_0)
# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end

# println("2")
# file_suffix = "sample_spring_one"
# p = [0.38 0.8]
# fps = 10
# Tf = 480
# destination = string("Nutrient Transport v.2/Simulation Results/Real/",file_suffix)



# global g; global avg_radii; global Len
# g, avg_radii, Len = import_MATLAB_network(file_suffix, true)
# params = graph_precomputation(g, avg_radii, Len)

# global Ne;  global Nv; global avg_radii; global Len; global V0; global kappa; global α; global D; global β; global Deff; global Dt; global τ; global C0; global E; global B; global K; global k; global L; global Ldag; global Mq_inv; global Minv1; global sum_Minv1; global elastic_modes; global elastic_eigvals; global edgelist; global out_neighbors; global in_neighbors; global Ein; global Eout; global EinT; global EoutT; global surface_area; global radius_sq; global Char_L; global Char_t; global R0; global Cq; global mu; global b 
# Ne, Nv, avg_radii, Len, V0, kappa, α, D, β, Deff, Dt, τ, C0, E, B, K, k, L, Ldag, Mq_inv, Minv1, sum_Minv1, elastic_modes, elastic_eigvals, edgelist, out_neighbors, in_neighbors, Ein, Eout, EinT, EoutT, surface_area, radius_sq, Char_L, Char_t, R0, Cq, mu, b = params


# phase_avg_disp = zeros(length(-pi:0.05:pi),1)
# global counter = 1
# for pt = -pi:0.05:pi
#     num_parts, num_per, dt, pos_0, ind_0, every_neighbor, Qin_t, Qout_t, A_t = modes_init_conds(Ne-1, Ne-2, pt, Ne, elastic_modes, fps, Tf)
#     pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, pos_0, ind_0)
#     avg_displacement = find_average_displacement(pos_time, ind_time)
#     global counter;
#     phase_avg_disp[counter] = maximum(avg_displacement)
#     counter = counter+1
#     println(counter)
# end

# output = string(destination,"phases2")
# open(output,"w") do io
#     writedlm(io,phase_avg_disp)
# end

# output = string(destination,"_avg_disp_mat2.txt")
# open(output,"w") do io
#     writedlm(io,avg_disp_mat)
# end


# avg_disp_mat = zeros(Ne-1,Ne-1)
# for i = 1:Ne-1
#     println(i)
#     for j = 1:Ne-1
#         if i >= j
#             global num_parts; global num_per; global dt; 
#             num_parts, num_per, dt, pos_0, ind_0, every_neighbor, Qin_t, Qout_t, A_t = modes_init_conds(Ne-i, Ne-j, 0.3687*pi, Ne, elastic_modes, fps, Tf)
#             pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, pos_0, ind_0)
#             avg_displacement = find_average_displacement(pos_time, ind_time)
#             avg_disp_mat[i,j] = maximum(avg_displacement)
#         end
#     end
# end


# test = zeros(60,1)
# for i = 1:60
  
#     global num_parts; global num_per; global dt; 
#     num_parts, num_per, dt, pos_0, ind_0, every_neighbor, Qin_t, Qout_t, A_t = modes_init_conds(Ne-1, Ne-i, 0.3687*pi, Ne, elastic_modes, fps, Tf)
#     pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, pos_0, ind_0)
#     avg_displacement = find_average_displacement(pos_time, ind_time)
#     test[i] = maximum(avg_displacement)
  
# end





# println("3")
# file_suffix = "sample_spring_three"
# p = [0.38 0.4]
# fps = 3
# sampleTime = 5
# destination = string("Nutrient Transport v.2/Simulation Results/Real/",file_suffix)

# @time num_parts, num_per, dt, pos_0, ind_0, every_neighbor, Qin_t, Qout_t, A_t = real_init_conds(file_suffix, fps, sampleTime)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, pos_0, ind_0)
# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end







# println("2")
# file_suffix = "sample_spring_eight"
# p = [0.38 0.8]
# fps = 5
# sampleTime = 3

# #model
# destination = string("Nutrient Transport v.2/Simulation Results/Model/",file_suffix)
# @time sol, elastic_modes, Ne = run_model(file_suffix, p, fps)

# @time A_t, Qin_t, Qout_t, pos_0, ind_0, num_parts, num_per, every_neighbor, dt, pos_0, ind_0 = model_init_conds(Ne, fps)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)
# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end

# destination = string("Nutrient Transport v.2/Simulation Results/Real/",file_suffix)
# #real
# @time num_parts, num_per, dt, pos_0, ind_0, every_neighbor, Qin_t, Qout_t, A_t = real_init_conds(file_suffix, fps, sampleTime)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)

# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end


# #plot
# @time plot_disp(file_suffix)




# println("3")
# file_suffix = "sample_spring_three"
# p = [0.38 0.4]
# fps = 3
# sampleTime = 5

# #model
# destination = string("Nutrient Transport v.2/Simulation Results/Model/",file_suffix)
# @time sol, elastic_modes, Ne = run_model(file_suffix, p, fps)

# @time A_t, Qin_t, Qout_t, pos_0, ind_0, num_parts, num_per, every_neighbor, dt, pos_0, ind_0 = model_init_conds(Ne, fps)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)
# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end

# destination = string("Nutrient Transport v.2/Simulation Results/Real/",file_suffix)
# #real
# @time num_parts, num_per, dt, pos_0, ind_0, every_neighbor, Qin_t, Qout_t, A_t = real_init_conds(file_suffix, fps, sampleTime)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)

# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end


# #plot
# @time plot_disp(file_suffix)




# println("4")
# file_suffix = "sample_spring_four"
# p = [0.38 0.3]
# fps = 3
# sampleTime = 5

# #model
# destination = string("Nutrient Transport v.2/Simulation Results/Model/",file_suffix)
# @time sol, elastic_modes, Ne = run_model(file_suffix, p, fps)

# @time A_t, Qin_t, Qout_t, pos_0, ind_0, num_parts, num_per, every_neighbor, dt, pos_0, ind_0 = model_init_conds(Ne, fps)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)
# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end

# destination = string("Nutrient Transport v.2/Simulation Results/Real/",file_suffix)
# #real
# @time num_parts, num_per, dt, pos_0, ind_0, every_neighbor, Qin_t, Qout_t, A_t = real_init_conds(file_suffix, fps, sampleTime)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)

# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end


# #plot
# @time plot_disp(file_suffix)




# println("5")
# file_suffix = "sample_spring_five"
# p = [0.38 0.7]
# fps = 3
# sampleTime = 5

# #model
# destination = string("Nutrient Transport v.2/Simulation Results/Model/",file_suffix)
# @time sol, elastic_modes, Ne = run_model(file_suffix, p, fps)

# @time A_t, Qin_t, Qout_t, pos_0, ind_0, num_parts, num_per, every_neighbor, dt, pos_0, ind_0 = model_init_conds(Ne, fps)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)
# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end

# destination = string("Nutrient Transport v.2/Simulation Results/Real/",file_suffix)
# #real
# @time num_parts, num_per, dt, pos_0, ind_0, every_neighbor, Qin_t, Qout_t, A_t = real_init_conds(file_suffix, fps, sampleTime)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)

# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end


# #plot
# @time plot_disp(file_suffix)




# println("6")
# file_suffix = "sample_spring_six"
# p = [0.38 0.6]
# fps = 3
# sampleTime = 5

# #model
# destination = string("Nutrient Transport v.2/Simulation Results/Model/",file_suffix)
# @time sol, elastic_modes, Ne = run_model(file_suffix, p, fps)

# @time A_t, Qin_t, Qout_t, pos_0, ind_0, num_parts, num_per, every_neighbor, dt, pos_0, ind_0 = model_init_conds(Ne, fps)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)
# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end

# destination = string("Nutrient Transport v.2/Simulation Results/Real/",file_suffix)
# #real
# @time num_parts, num_per, dt, pos_0, ind_0, every_neighbor, Qin_t, Qout_t, A_t = real_init_conds(file_suffix, fps, sampleTime)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)

# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end


# #plot
# @time plot_disp(file_suffix)




# println("7")
# file_suffix = "sample_spring_seven"
# p = [0.38 0.4]
# fps = 3
# sampleTime = 5

# #model
# destination = string("Nutrient Transport v.2/Simulation Results/Model/",file_suffix)
# @time sol, elastic_modes, Ne = run_model(file_suffix, p, fps)

# @time A_t, Qin_t, Qout_t, pos_0, ind_0, num_parts, num_per, every_neighbor, dt, pos_0, ind_0 = model_init_conds(Ne, fps)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)
# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end

# destination = string("Nutrient Transport v.2/Simulation Results/Real/",file_suffix)
# #real
# @time num_parts, num_per, dt, pos_0, ind_0, every_neighbor, Qin_t, Qout_t, A_t = real_init_conds(file_suffix, fps, sampleTime)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)

# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end


# #plot
# @time plot_disp(file_suffix)




# println("8")
# file_suffix = "sample_spring_eight"
# p = [0.38 0.2]
# fps = 3
# sampleTime = 5

# #model
# destination = string("Nutrient Transport v.2/Simulation Results/Model/",file_suffix)
# @time sol, elastic_modes, Ne = run_model(file_suffix, p, fps)

# @time A_t, Qin_t, Qout_t, pos_0, ind_0, num_parts, num_per, every_neighbor, dt, pos_0, ind_0 = model_init_conds(Ne, fps)
# pos_0[1:Ne] = 0.5*Len
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)
# # @time avg_displacement = find_average_displacement(pos_time, ind_time)

# # output = string(destination,"_avg_disp.txt")
# # open(output,"w") do io
# #     writedlm(io,avg_displacement)
# # end

# ###  SAVING particle positions

# output = string("init_part_pos.txt")
# open(output,"w") do io
#     writedlm(io,pos_0[1:Ne] ./ Len)
# end

# fin_part_pos = zeros(Ne,2)
# for i = 1:Ne
#     fin_part_pos[i,1] = pos_time[end,i] / Len[ind_time[end,i]]
#     fin_part_pos[i,2] = ind_time[end,i]
# end

# output = string("fin_part_pos.txt")
# open(output, "w") do io
#     writedlm(io,fin_part_pos)
# end


# #convert to amp.
# amp_0 = (sol.u[1][1:Ne] .- V0) ./ V0
# amp_f = (sol.u[end][1:Ne] .- V0) ./ V0

# output = string("amp0.txt")
# open(output, "w") do io
#     writedlm(io,amp_0)
# end

# output = string("ampf.txt")
# open(output, "w") do io
#     writedlm(io,amp_f)
# end


# destination = string("Nutrient Transport v.2/Simulation Results/Real/",file_suffix)
# #real
# @time num_parts, num_per, dt, pos_0, ind_0, every_neighbor, Qin_t, Qout_t, A_t = real_init_conds(file_suffix, fps, sampleTime)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)

# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end


# #plot
# @time plot_disp(file_suffix)





# println("9")
# file_suffix = "sample_spring_nine"
# p = [0.38 0.3]
# fps = 3
# sampleTime = 5

# #model
# destination = string("Nutrient Transport v.2/Simulation Results/Model/",file_suffix)
# @time sol, elastic_modes, Ne = run_model(file_suffix, p, fps)

# @time A_t, Qin_t, Qout_t, pos_0, ind_0, num_parts, num_per, every_neighbor, dt, pos_0, ind_0 = model_init_conds(Ne, fps)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)
# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end

# destination = string("Nutrient Transport v.2/Simulation Results/Real/",file_suffix)
# #real
# @time num_parts, num_per, dt, pos_0, ind_0, every_neighbor, Qin_t, Qout_t, A_t = real_init_conds(file_suffix, fps, sampleTime)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)

# @time avg_displacement = find_average_displacement(pos_time, ind_time)

# output = string(destination,"_avg_disp.txt")
# open(output,"w") do io
#     writedlm(io,avg_displacement)
# end


# #plot
# @time plot_disp(file_suffix)


# fps = 10
# Tf = 75
# num_parts = 3Ne

# global pos_0;
# pos_0 = []
#     for i = 1:num_per
#         append!(pos_0,i.*Len/(num_per+1))
#     end

#     global ind_0;
#     ind_0 = []
#     for i = 1:num_per
#         append!(ind_0,[i for i = 1:Ne])
#     end


# @time Qin_t, Qout_t, A_t = modes_init_conds(Ne-2, Ne-1, 0, fps, Tf)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)
# @time disp0 = find_average_displacement(pos_time, ind_time)

# @time Qin_t, Qout_t, A_t = modes_init_conds(Ne-2, Ne-1, 1pi/6, fps, Tf)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)
# @time disp1 = find_average_displacement(pos_time, ind_time)

# @time Qin_t, Qout_t, A_t = modes_init_conds(Ne-2, Ne-1, 2pi/6, fps, Tf)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)
# @time disp2 = find_average_displacement(pos_time, ind_time)

# @time Qin_t, Qout_t, A_t = modes_init_conds(Ne-2, Ne-1, 3pi/6, fps, Tf)
# @time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, num_parts)
# @time disp3 = find_average_displacement(pos_time, ind_time)


println("8")
file_suffix = "sample_spring_eight"
#@time sol, elastic_modes, Ne = run_model(file_suffix, p, fps)
p = [0.38 0.2]
fps = 3
sampleTime = 5
fps = 10
Tf = 1000
dt = 2/fps

@time Qin_t, Qout_t, A_t, radii_time_modes, pos_0, ind_0, every_neighbor = modes_init_conds(Ne-2, Ne-1, pi/2, fps, Tf)
@time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, 3Ne)


output = string("init_part_pos.txt")
open(output,"w") do io
    writedlm(io,pos_0[1:Ne] ./ Len)
end

fin_part_pos = zeros(Ne,2)
for i = 1:Ne
    fin_part_pos[i,1] = pos_time[end,i] / Len[ind_time[end,i]]
    fin_part_pos[i,2] = ind_time[end,i]
end

output = string("fin_part_pos.txt")
open(output, "w") do io
    writedlm(io,fin_part_pos)
end


#convert to amp.
for j = 1:Ne
    temp = radii_time_modes[:,j] 
    radii_time_modes[:,j] = (temp .- mean(temp)) ./ (maximum(temp) - mean(temp))
end

amp_0 = radii_time_modes[1,:]
amp_f = radii_time_modes[end-10,:]

output = string("amp0.txt")
open(output, "w") do io
    writedlm(io,amp_0)
end

output = string("ampf.txt")
open(output, "w") do io
    writedlm(io,amp_f)
end
