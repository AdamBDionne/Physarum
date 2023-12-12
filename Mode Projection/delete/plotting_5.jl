using StatsBase
using LightGraphs 
using DifferentialEquations
using Plots
using Images



function project_real_and_model2(file_suffix, sol, elastic_modes, Ne, fps)

    #put model signal in same format 
    radii_time_model = zeros(Ne, size(sol.u,1))
    numSamples = size(radii_time_model,2)
    for j = 1:numSamples
        radii_time_model[:,j] = sol.u[j][1:Ne]
    end
    for i = 1:Ne
        radii_time_model[i,:] = (radii_time_model[i,:] .- mean(radii_time_model[i,:])) ./ maximum(abs.(radii_time_model[i,:] .- mean(radii_time_model[i,:])))
    end
    
    modes = copy(elastic_modes)

    i = Ne-1
    time_series_model = zeros(numSamples)
    #projection time series
    for j = 1:numSamples
        norm_model = ( (sum(radii_time_model[:,j].^2)) * (sum(modes[:,i].^2)) )^(1/2)
        time_series_model[j] = sum(radii_time_model[:,j] .* modes[:,i]) / (norm_model)
    end


    curr_max = maximum(time_series_model)
    frame_model = radii_time_model[:,argmax(time_series_model)]
    return curr_max, frame_model
end



file_suffix = "sample_spring_seven"
p = [0.38 0.4]
fps = 3

global u0
u0 = readdlm("test_ode.txt")
u0 = u0[end,:]

# file_suffix = "sample_spring_six"
# p = [0.38 0.6]
global cmax = 0.0
global u0M
save = false 
maxes = []
for tst = 1:1
    global u0
    global cmax; global frame_real; global frame_model; global frame_conc; global save; global u0M
   
    sol, elastic_modes, Ne, u0 = @time run_model(file_suffix, p, fps, u0)
    curr_max, c_frame_model = project_real_and_model2(file_suffix, sol, elastic_modes, Ne, fps)
    push!(maxes, curr_max)
    println(curr_max)
    if curr_max > cmax
        frame_model = c_frame_model;
        cmax = curr_max 
        println("New MAX")
        save = true;
        u0M = u0
        println(cmax)
    end
end

save = false

if save

    output = "frame_model2.txt"
    open(output,"w") do io
        writedlm(io,frame_model)
    end

    output = "frame_model2_u0.txt"
    open(output, "w") do io
        writedlm(io, u0M)
    end

end
