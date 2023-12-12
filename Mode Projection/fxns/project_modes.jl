## Project real and model amplitudes onto modes
function project_real_and_model(file_suffix, sol, elastic_modes, Ne, fps)
    source = string("Data/",file_suffix,"/",file_suffix)
    #destination = string("Processed data/",file_suffix,"/",file_suffix)

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
        #radii_time[j,:] = (radii_data_curr .- mean(radii_data_curr)) ./ mean(radii_data_curr)
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
        #radii_time_model[i,:] = sqrt.(radii_time_model ./ (pi * Len[i]))
        radii_time_model[i,:] = (radii_time_model[i,:] .- mean(radii_time_model[i,:])) ./ maximum(abs.(radii_time_model[i,:] .- mean(radii_time_model[i,:])))
        #radii_time_model[i,:] = (radii_time_model[i,:] .- mean(radii_time_model[i,:])) ./ mean(radii_time_model[i,:])
    end

    #project both
    numSamples = size(radii_time,2)
    avg_amp = zeros(Ne-1,1)
    time_series = zeros(Ne-1,numSamples)
    time_series_rms = zeros(Ne-1,numSamples)

    avg_amp_model = zeros(Ne-1,1)
    time_series_model = zeros(Ne-1,numSamples)
    time_series_rms_model = zeros(Ne-1,numSamples)

    for i = 1:Ne-1
        #projection time series
        for j = 1:numSamples
            norm = ( (sum(radii_time[:,j].^2)) * (sum(modes[:,i].^2)) )^(1/2)
            norm_model = ( (sum(radii_time_model[:,j].^2)) * (sum(modes[:,i].^2)) )^(1/2)
            time_series[i,j] = sum(radii_time[:,j] .* modes[:,i]) ./ (norm)
            time_series_model[i,j] = sum(radii_time_model[:,j] .* modes[:,i]) / (norm_model)
        end


        #rms average time_series
        for j = 1:numSamples
            time_series_rms[i,j] = sqrt.(mean((time_series[i,maximum([1;j-45]):minimum([numSamples;j+45])]).^2))
            time_series_rms_model[i,j] = sqrt.(mean((time_series_model[i,maximum([1;j-45]):minimum([numSamples;j+45])]).^2))
        end

        for j = 1:numSamples
            time_series_rms[i,j] = sqrt.(mean((time_series_rms[i,maximum([1;j-45]):minimum([numSamples;j+45])]).^2))
            time_series_rms_model[i,j] = sqrt.(mean((time_series_rms_model[i,maximum([1;j-45]):minimum([numSamples;j+45])]).^2))
        end
       
        avg_amp[i] = maximum(time_series_rms[i,:])
        avg_amp_model[i] = maximum(time_series_rms_model[i,:])
    end

    # #save modes
    # output = string(destination,"_MODES.txt")
    # open(output,"w") do io
    #     writedlm(io,modes)
    # end

        # #print two maximum
        # max = argmax(avg_amp)
        # max = max[1]
        # temp = copy(avg_amp); temp[max] = 0;
        # second_max = argmax(temp)
        # second_max = second_max[1]

        # #find phase difference between two max for real and model
        # sample_1 = time_series[max,:]
        # sample_2 = time_series[second_max,:]
        # #phase_dif = acos(sum(sample_1 .* sample_2) / sqrt(sum(sample_1.^2) * sum(sample_2.^2)) ) 

        # #find time period 
        # t = crosscor(sample_1, sample_1, [i for i = 1:200])
        # maxima = findlocalmaxima(t)
        # maxima = [getindex(max, 1) for max in maxima]
        # tp_sample_1 = maxima[2] / fps

        # #phase difference
        # tx = crosscor(sample_1, sample_2, [i for i = 1:200])
        # maxima = findlocalmaxima(tx)
        # maxima = [getindex(max, 1) for max in maxima]
        # if maxima[1] == 1
        #     t_diff = (maxima[2]) / fps
        # else 
        #     t_diff = (maxima[1]) / fps
        # end
        # phase_dif = 2pi - ((t_diff / tp_sample_1) * 2pi) 

        # sample_1_model = time_series_model[max,:]
        # sample_2_model = time_series_model[second_max,:]
        # #phase_dif_model = acos(sum(sample_1_model .* sample_2_model) / sqrt(sum(sample_1_model.^2) * sum(sample_2_model.^2)) ) 
        # #find time period 
        # t = crosscor(sample_1_model, sample_1_model, [i for i = 1:200])
        # maxima = findlocalmaxima(t)
        # maxima = [getindex(max, 1) for max in maxima]
        # tp_sample_1 = maxima[2] / fps

        # #phase difference
        # tx = crosscor(sample_1_model, sample_2_model, [i for i = 1:200])
        # maxima = findlocalmaxima(tx)
        # maxima = [getindex(max, 1) for max in maxima]
        # if maxima[1] == 1
        #     t_diff = (maxima[2]) / fps
        # else 
        #     t_diff = (maxima[1]) / fps
        # end
        # phase_dif_model = 2pi - ((t_diff / tp_sample_1) * 2pi) 


        # output = string(destination,"_mode_proj_info.txt")
        # open(output,"w") do io
        #     writedlm(io,[max;second_max;phase_dif;phase_dif_model])
        # end

    #find maximum proj for real 

    indmaxr = argmax(avg_amp)[1]
    indmaxm = argmax(avg_amp_model)[1]
    ind_real = argmax(time_series_rms)[2]
    ind_model = argmax(time_series_rms_model)[2]
    
  #  ind_real = argmax(time_series[ind_max_real,:])[1]
  #  ind_model = argmax(time_series_model[ind_max_real,:])[1]

    #save the frames for each

        # #save the times
        # output = "time.txt"
        # open(output,"w") do io
        #     writedlm(io,[3*ind_real/60, 3*ind_model/60])
        # end

        # output = "frame_real.txt"
        # open(output,"w") do io
        #     writedlm(io,radii_time[:,ind_real])
        # end

    output = "Mode projection/results/model_largest_proj.txt"
    open(output,"w") do io
        writedlm(io,radii_time_model[:,ind_model])
    end

        # output = "model_conc.txt"
        # open(output,"w") do io
        #     writedlm(io,sol.u[ind_real][Ne+1:end])
        # end


    output = "Mode projection/results/mode_largest_proj.txt"
    open(output,"w") do io
        writedlm(io,modes[:,indmaxm])
    end

    #return abs.(time_series), abs.(time_series_model), avg_amp, avg_amp_model, ind_real, ind_model

    return time_series_rms, time_series_rms_model, avg_amp, avg_amp_model, ind_real, ind_model, indmaxm, indmaxr
end

