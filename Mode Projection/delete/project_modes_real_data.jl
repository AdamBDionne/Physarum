#Project modes onto model solution
using Images
using PyPlot
using PyCall

file_suffix = "sample_spring_three"
source = string("Raw data/",file_suffix,"/",file_suffix)

#include("Phase Diagrams/fxns/compute_jacobian_shifted.jl")
g, avg_radii, Len = import_MATLAB_network(file_suffix, true)
params = graph_precomputation(g, avg_radii, Len)
Ne = ne(g)
translator = readdlm(string(source,"_translator.txt"))
translator = Int.(abs.(translator))

filename_radii = string(source,"_radii_time.txt")
radii_data = readdlm(filename_radii) #.* 1.83 * 10^(-6)
#find number of samples 
broken_up = split(radii_data[1],",")
radii_time = zeros(ne(g),size(broken_up,1))

for j = 1:ne(g)
    broken_up = split(radii_data[j],",")
    radii_data_curr = parse.(Float64,broken_up)
    radii_time[j,:] = radii_data_curr
end
translator2 = 0*translator;
for i =1:Ne
    translator2[translator[i]] = i 
end
radii_time = radii_time[translator[:,1],:]

# Save parameters
Ne, Nv, avg_radii, Len, V0, kappa, α, D, β, Deff, Dt, τ, C0, E, B, K, k, L, Ldag, Mq_inv, Minv1, sum_Minv1, elastic_modes, elastic_eigvals, edgelist, out_neighbors, in_neighbors, Ein, Eout, EinT, EoutT, ϵ_s, ϵ_c, surface_area, radius_sq, Char_L = params

#Normalize modes
modes = copy(elastic_modes)
for j = 1:Ne
    modes[:,j] = elastic_modes[:,j] .- mean(elastic_modes[:,j])
    modes[:,j] = modes[:,j] ./ maximum(abs.(modes[:,j]))  

    radii_time[j,:] = radii_time[j,:] .- mean(radii_time[j,:])
    radii_time[j,:] = radii_time[j,:] ./ maximum(abs.(radii_time[j,:]))
end

#save modes (translated)
#modes_translated = modes[:,translator2[:,1]]
output = string(source, "_MODES.txt")
open(output,"w") do io
    writedlm(io,modes)
end


numSamples = size(radii_time,2)
avg_amp = zeros(size(modes)[1]-1,1)
time_series = zeros(Ne,numSamples)
time_series_rms = zeros(Ne,numSamples)

for i = 1:size(modes)[1]-1
    #projection time series
    for j = 1:numSamples
        norm = ( (sum(radii_time[:,j].^2)) * (sum(modes[:,i].^2)) )^(1/2)
        time_series[i,j] = sum(radii_time[:,j] .* modes[:,i]) ./ (norm)
    end

    #average amplitude in time series
    avg_amp[i] = mean(abs.(time_series[i,findlocalmaxima(abs.(time_series[i,:]))]))

    #rms average time_series
    for j = 1:numSamples
        time_series_rms[i,j] = mean(abs.(time_series[i,maximum([1;j-300]):minimum([numSamples;j+300])]))
    end
end

sums_ordered = zeros(Ne-1,2)
for i = 1:Ne-1
    sums_ordered[i,1] = i;
    sums_ordered[i,2] = avg_amp[i];
end

sums_ordered = sums_ordered[sortperm(sums_ordered[:,2]),:]

#order parameter

#get phase vs time 
phase_time = zeros(Ne, numSamples)
for i = 1:Ne
    phase_time[i,:] = angle.(hilbert(reshape(radii_time[i,:],numSamples,1)))
end

global order_time_series = zeros(numSamples,1)
global counter = 0;
for i = 1:Nv
    out = out_neighbors[i]
    inn = in_neighbors[i]
    neighbs = [out ; inn]
    adjs = collect(combinations(neighbs,2))
    for j = 1:length(adjs)
        phase_dif = angle.(phase_time[adjs[j][1],:] ./ phase_time[adjs[j][2],:])
        global order_time_series = order_time_series .+ exp.(im*(phase_dif))
        global counter = counter + 1;
    end
end
order_time_series = real.(order_time_series) ./ counter

#phase difference between two dominant projections
sample_1 = time_series[Int(sums_ordered[end,1]),:]
sample_2 = time_series[Int(sums_ordered[end-1,1]),:]

    # #hilberts
    # sample_1_h = hilbert(reshape(sample_1,numSamples,1))
    # sample_2_h = hilbert(reshape(sample_2,numSamples,1))

    # phase_dif = angle.(sample_1_h ./ sample_2_h)

phase_dif = acos(sum(sample_1 .* sample_2) / sqrt(sum(sample_1.^2) * sum(sample_2.^2)) ) 



# # Plot results 

# PyPlot.close()
# rcParams = PyDict(matplotlib["rcParams"])
# rcParams["font.size"] = 18

# fig, ax = PyPlot.subplots(figsize=(18,9))

# time = [(1/60)*i*10/fps for i in t_i:1/fps:t_f]

# #plot time series
# plt.subplot(1,4,1)
# for i = 1:size(modes)[1]-1
#     PyPlot.plot(time, abs.(time_series[i,:]))
# end
# PyPlot.xlabel("Time (seconds)")
# PyPlot.ylabel("Mode projection")

# #plot rms
# plt.subplot(1,4,2)
# for i = 1:size(modes)[1]-1
#     PyPlot.plot(time, time_series_rms[i,:])
# end
# PyPlot.xlabel("Time (seconds)")
# PyPlot.ylabel("Mode projection RMS")

# #plot order of proj
# plt.subplot(1,4,3)
# PyPlot.scatter([i for i = 1:size(modes)[1]-1], avg_amp)
# PyPlot.xlabel("Mode Number")
# PyPlot.ylabel("Average Projection Amplitude")

# #plot two largest time series
# plt.subplot(1,4,4)
# PyPlot.plot(time, time_series[Int(sums_ordered[end,1]),:])
# PyPlot.plot(time, time_series[Int(sums_ordered[end-1,1]),:])
# PyPlot.xlabel("Time (seconds)")
# PyPlot.ylabel("Mode Projection")


# plt.tight_layout()
# PyPlot.savefig("mode_projections.png")