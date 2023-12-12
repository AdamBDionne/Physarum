#Project modes onto model solution
using Images
using PyPlot
using PyCall

modes = copy(elastic_modes)
for j = 1:Ne
    modes[:,j] = elastic_modes[:,j] .- mean(elastic_modes[:,j])
    modes[:,j] = modes[:,j] ./ maximum(abs.(modes[:,j]))  
end

numSamples = length(sol.u)
avg_amp = zeros(size(modes)[1]-1,1)
time_series = zeros(Ne,numSamples)
time_series_rms = zeros(Ne,numSamples)

for i = 1:size(modes)[1]-1
    #projection time series
    for j = 1:numSamples
        norm = ( (sum(sol.u[j][1:Ne].^2)) * (sum(modes[:,i].^2)) )^(1/2)
        time_series[i,j] = sum(sol.u[j][1:Ne] .* modes[:,i]) ./ (norm)
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


# Plot results 

PyPlot.close()
rcParams = PyDict(matplotlib["rcParams"])
rcParams["font.size"] = 18

fig, ax = PyPlot.subplots(figsize=(18,9))

time = [(1/60)*i*10/fps for i in t_i:1/fps:t_f]

#plot time series
plt.subplot(1,4,1)
for i = 1:size(modes)[1]-1
    PyPlot.plot(time, abs.(time_series[i,:]))
end
PyPlot.xlabel("Time (seconds)")
PyPlot.ylabel("Mode projection")

#plot rms
plt.subplot(1,4,2)
for i = 1:size(modes)[1]-1
    PyPlot.plot(time, time_series_rms[i,:])
end
PyPlot.xlabel("Time (seconds)")
PyPlot.ylabel("Mode projection RMS")

#plot order of proj
plt.subplot(1,4,3)
PyPlot.scatter([i for i = 1:size(modes)[1]-1], avg_amp)
PyPlot.xlabel("Mode Number")
PyPlot.ylabel("Average Projection Amplitude")

#plot two largest time series
plt.subplot(1,4,4)
PyPlot.plot(time, time_series[Int(sums_ordered[end,1]),:])
PyPlot.plot(time, time_series[Int(sums_ordered[end-1,1]),:])
PyPlot.xlabel("Time (seconds)")
PyPlot.ylabel("Mode Projection")


plt.tight_layout()
PyPlot.savefig("mode_projections.png")