## Simulate nutrient transport on mode super position at phi = pi/2 phase shift
include("fxns/helper_functions.jl")
include("fxns/model_helper_functions.jl")
include("fxns/simulate_modes.jl")
using MAT


file_suffix = "sample two"
destination = "Nutrient Transport/results/"

global g; global avg_radii; global Len
g, avg_radii, Len = import_MATLAB_network(file_suffix, true)
params = graph_precomputation(g, avg_radii, Len)

global Ne;  global Nv; global avg_radii; global Len; global V0; global kappa; global α; global D; global β; global Deff; global Dt; global τ; global C0; global E; global B; global K; global k; global L; global Ldag; global Mq_inv; global Minv1; global sum_Minv1; global elastic_modes; global elastic_eigvals; global edgelist; global out_neighbors; global in_neighbors; global Ein; global Eout; global EinT; global EoutT; global surface_area; global radius_sq; global Char_L; global Char_t; global R0; global Cq; global mu; global b 
Ne, Nv, avg_radii, Len, V0, kappa, α, D, β, Deff, Dt, τ, C0, E, B, K, k, L, Ldag, Mq_inv, Minv1, sum_Minv1, elastic_modes, elastic_eigvals, edgelist, out_neighbors, in_neighbors, Ein, Eout, EinT, EoutT, surface_area, radius_sq, Char_L, Char_t, R0, Cq, mu, b = params

sampleTime = 5
fps = 10
Tf = 1000
dt = 2/fps

@time Qin_t, Qout_t, A_t, radii_time_modes, pos_0, ind_0, every_neighbor = modes_init_conds(Ne-2, Ne-1, pi/2, fps, Tf)
@time pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, 3Ne)


#convert to amp.
for j = 1:Ne
    temp = radii_time_modes[:,j] 
    radii_time_modes[:,j] = (temp .- mean(temp)) ./ (maximum(temp) - mean(temp))
end

## saving results for movie 
num_frames = 900
step_size = floor(Int, size(radii_time_modes,1) / num_frames)
step_size += step_size % 2 
start_index = 2  
indices1 = start_index:step_size:size(radii_time_modes,1)
indices2 = 1:round(Int,step_size/2):round(Int,size(radii_time_modes,1)/2)
indices1 = 1:2:1800
indices2 = 1:900

output = string(destination, "amp_t.txt")
open(output, "w") do io
    writedlm(io,radii_time_modes[indices1,:])
end

fin_part_pos = zeros(Ne,2,length(indices2))
c = 1
for k = indices2
    for i = 1:Ne
        fin_part_pos[i,1,c] = pos_time[k,i] / Len[ind_time[k,i]]
        fin_part_pos[i,2,c] = ind_time[k,i]
    end
    c = c+1
end

file = matopen(string(destination, "fin_part_pos.mat"), "w")
write(file, "fin_part_pos", fin_part_pos)
close(file)
