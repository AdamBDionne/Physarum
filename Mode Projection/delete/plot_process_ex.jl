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
        #radii_time[j,:] = (radii_data_curr .- mean(radii_data_curr)) ./ maximum(abs.(radii_data_curr .- mean(radii_data_curr)))
        radii_time[j,:] = (radii_data_curr .- mean(radii_data_curr)) ./ mean(radii_data_curr)
    end


    #translate from MATLAB listing to Julia listing
    translator2 = 0*translator;
    for i =1:Ne
        translator2[translator[i]] = i 
    end
    radii_time = radii_time[translator[:,1],:]

    #Normalize modes
    modes = copy(elastic_modes)

    #put model signal in same format 
    radii_time_model = 0 .* radii_time
    for j = 1:size(radii_time,2)
        radii_time_model[:,j] = sol.u[j][1:Ne]
    end
    for i = 1:Ne
        radii_time_model[i,:] = sqrt.(radii_time_model[i,:] ./ (pi * Len[i]))
        #radii_time_model[i,:] = (radii_time_model[i,:] .- mean(radii_time_model[i,:])) ./ maximum(abs.(radii_time_model[i,:] .- mean(radii_time_model[i,:])))
        radii_time_model[i,:] = (radii_time_model[i,:] .- mean(radii_time_model[i,:])) ./ mean(radii_time_model[i,:])
    end

    #project on eachother
    numSamples = size(radii_time,2)
    proj_together = zeros(numSamples)

    for i = 1:Ne-1
        #projection time series
        for j = 1:numSamples
            norm = ( (sum(radii_time[:,j].^2)) * (sum(radii_time_model[:,j].^2)) )^(1/2)
            proj_together[j] = sum(radii_time[:,j] .* radii_time_model[:,j]) ./ (norm)
        end
    end

    frame_max = argmax(proj_together)
    curr_max = proj_together[frame_max]
    c_frame_real = radii_time[:,frame_max] 
    c_frame_model = radii_time[:,frame_max]
    c_frame_conc = sol.u[frame_max][Ne+1:end]
    
    return curr_max, c_frame_real, c_frame_model, c_frame_conc
end

