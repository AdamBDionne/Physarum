using StatsBase
using LightGraphs 
using DifferentialEquations
using Plots
using Images



function project_real_and_model2(file_suffix, sol, elastic_modes, Ne, fps)

    source = string("Raw data/",file_suffix,"/",file_suffix)
    destination = string("Processed data/",file_suffix,"/",file_suffix)

    translator = readdlm(string(source,"_translator.txt"))
    translator = Int.(abs.(translator))

    filename_radii = string(source,"_radii_time.txt")
    radii_data = readdlm(filename_radii)

    #find number of samples 
    broken_up = split(radii_data[1],",")
    radii_time = zeros(Ne,size(broken_up,1))

    #prase and normalize
    for j = 1:Ne
        broken_up = split(radii_data[j],",")
        radii_data_curr = parse.(Float64,broken_up)
        radii_time[j,:] = (radii_data_curr .- mean(radii_data_curr)) ./ maximum(abs.(radii_data_curr .- mean(radii_data_curr)))
    end

    #translate from MATLAB listing to Julia listing
    translator2 = 0*translator;
    for i =1:Ne
        translator2[translator[i]] = i 
    end
    radii_time = radii_time[translator[:,1],:]

    #put model signal in same format 
    radii_time_model = 0 .* radii_time
    for j = 1:size(radii_time,2)
        radii_time_model[:,j] = sol.u[j][1:Ne]
    end
    for i = 1:Ne
        radii_time_model[i,:] = (radii_time_model[i,:] .- mean(radii_time_model[i,:])) ./ maximum(abs.(radii_time_model[i,:] .- mean(radii_time_model[i,:])))
    end

    norm = abs.(radii_time') * abs.(radii_time_model)
    projects = radii_time' * radii_time_model ./ norm
    max_proj_ind = argmax(projects)

    ind_real = max_proj_ind[1]
    ind_model = max_proj_ind[2]

    curr_max = maximum(projects)
    frame_real = radii_time[:,ind_real]
    frame_model = radii_time_model[:,ind_model]

    t = permutedims(hcat(sol.u...))
    frame_conc = sol.u[ind_model][Ne+1:end] ./ mean(t,dims=1)[Ne+1:end]

    return curr_max, frame_real, frame_model, frame_conc
end



file_suffix = "sample_spring_three"
p = [0.38 0.4]
fps = 3

# file_suffix = "sample_spring_six"
# p = [0.38 0.6]
global cmax = 0.2
save = true
maxes = []
for tst = 1:10
    global cmax; global frame_real; global frame_model; global frame_conc; global save
    println(tst)
    sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
    curr_max, c_frame_real, c_frame_model, c_frame_conc = project_real_and_model2(file_suffix, sol, elastic_modes, Ne, fps)
    push!(maxes, curr_max)
 
    if curr_max > cmax
        frame_real = c_frame_real; frame_model = c_frame_model;
        frame_conc = c_frame_conc; cmax = curr_max 
        println("New MAX")
        save = true;
        println(cmax)
    end
end

if save
    output = "frame_real_a.txt"
    open(output,"w") do io
        writedlm(io,frame_real)
    end

    output = "frame_model_a.txt"
    open(output,"w") do io
        writedlm(io,frame_model)
    end

    output = "model_conc_a.txt"
    open(output,"w") do io
        writedlm(io,frame_conc)
    end
end
