using Statistics
using FFTW
using StatsBase
using Images

#find time period, phase difference between first and second harmonic,
#and amplitude
function fft_info(sample,sampleTime)
    time_periods = zeros(Ne)
    amplitudes = zeros(Ne)
    THDF = zeros(Ne)

    #go through each edge and find time period and amplitude
    #we will then average these 

    for edge_ind = 1:Ne
        #get signal
        signal = zeros(Float64,size(sample)[1])
        for i in 1:size(sample)[1]
            signal[i] = sample[i][edge_ind]
        end
        avg = mean(signal)
        rad_sig = sqrt.(signal ./ (pi * Len[edge_ind]))
        signal = signal .- mean(signal)

        #cross correlation 
        t = crosscor(signal, signal, [i for i = 1:1000])
        maxima = findlocalmaxima(t)
        maxima = [getindex(max, 1) for max in maxima]
        time_periods[edge_ind] = maxima[2] * sampleTime

        #find signal's local maximums 
        maxima = findlocalmaxima(signal)
        maxima = [getindex(max, 1) for max in maxima]
        num_maxima = length(maxima)
        if maxima[1] > 10
            ind_start = maxima[1] 
            max_start = 1
        else
            ind_start = maxima[2]
            max_start = 2
        end
        if maxima[end] < num_maxima-10
            ind_fin = maxima[end]
            max_fin = num_maxima
        else
            ind_fin = maxima[end-1]
            max_fin = num_maxima-1
        end

        #use these to find amplitude 
        
        amplitudes[edge_ind] = median(rad_sig[maxima[max_start:max_fin]]) / mean(rad_sig)

        signal_max_to_max = signal[ind_start:ind_fin]
        N = length(signal_max_to_max)

        F = fft(signal_max_to_max) |> fftshift
        freqs = fftfreq(N, 1.0/sampleTime) |> fftshift

        zero_ind = argmax(freqs .== 0)
        F_pos = F[zero_ind:end]
        freqs_pos = freqs[zero_ind:end]

        #find maxima 
        maxima = partialsortperm(abs.(F_pos[5:60]),1:3,rev=true)

        first_harmonic_freq = freqs_pos[maxima[1]]

        #time_periods[edge_ind] = 1/first_harmonic_freq

        #quantify how nonlinear this signal is
            #num_harmonics = floor(size(freqs_pos,1) / maxima[1])
            #THDF = sqrt(sum([abs(F_pos[i*maxima[1]])^2 for i=2:Int(num_harmonics)])) / abs(F_pos[maxima[1]])
    end

    return median(time_periods), median(amplitudes)
end

#find the local and global order of sample 
function find_order(sample)

    num_samples = size(sample,1)
    phases = zeros(Ne, num_samples)

    for index in 1:Ne
        #get vessel's signal 
        signal = zeros(Float64,num_samples)
        for i in 1:size(sample)[1]
            signal[i] = sample[i][index]
        end
        signal = signal .- mean(signal)

        #find analytic signal
        signal_hil = hilbert(reshape(signal,num_samples,1))
        phases[index,:] = angle.(signal_hil)
    end

    #global order
    global_order_signal = zeros(num_samples,1)
    for i = 1:Ne
        global_order_signal += exp.(im .* phases[i,:])
    end

    global_order_signal = abs.(global_order_signal)
    global_order_signal = global_order_signal ./ Ne

    #local order
    local_order_signal = zeros(num_samples,1)
    for i = 1:Nv
        out = out_neighbors[i]
        inn = in_neighbors[i]
        neighbs = [out ; inn]
        curr_order = [0+0im for j = 1:num_samples]
        for e in neighbs
            curr_order += exp.(im .* phases[e,:])
        end
        local_order_signal += abs.(curr_order) ./ length(neighbs)
    end
    local_order_signal = local_order_signal ./ Nv

    return median(global_order_signal), median(local_order_signal)
end

function find_flow_amp(sample, sampleTime, p)
    #let's sample over 20 simulation seconds 
    sample_inds = [i for i = 1:Int(20/sampleTime)]
    edge_max_flow_vel = zeros(Ne)
    for ind in sample_inds
        flow_vels = find_flow_vels(sample[ind],p)
        for j in 1:Ne
            curr_flow_vel = abs(flow_vels[j])
            if curr_flow_vel > edge_max_flow_vel[j]
                edge_max_flow_vel[j] = curr_flow_vel
            end
        end
    end
    median_vel = median(edge_max_flow_vel)
    #let's convert to um / s 
    unitful_vel = ( median_vel * (Char_L/Char_t) ) * 10^6
    return 2*unitful_vel, std(edge_max_flow_vel .* (Char_L/Char_t) .* 10^6)
end

#find flows
function find_flow_vels(u, p)
    ϵ_s, ϵ_c = p
    V = u[1:Ne] #volume 
    C = u[Ne+1:2Ne] #solute amount
    radius_sq = (V) ./ (pi .* Len)
    R = sqrt.(radius_sq)

    ϵ = (R .- R0)./R0 #percent change from equilibrium volume

    
    # nonlinear conductivities and damping:
        # K = (pi.*R.^4)./(8*mu.*Len)
        # τ = diagm(0 => vec((b ./(4 .* pi^2 .* radius_sq .* Len.^2)) .+ (1 ./ (12*K))))
        # Mq = τ + Cq
        # Mq_inv = factorize(Mq)
        # Minv1 = Mq_inv \ ones(Ne)
        # sum_Minv1 = sum(Minv1)

    #Volume dynamics
                  #  restoring force      nonlinear fudge    actomyosin contractions 
    r = Mq_inv \ (-(k .*ϵ) .- (kappa .*(ϵ).^3) .- ( α .*(C./C0).*(1 .- ϵ/ϵ_s)))
    μ = sum(r) / sum_Minv1
    Vdot = r .- (Minv1 .* μ)

    #Flows 
    volFlows = mid_flows(Vdot)
    return volFlows ./ (pi * R.^2)
end


function save_phase_info(result,filename,curr_size)
    output = string("Phase Diagrams/",filename,"/",curr_size)
    for j = 1:14
        savefig(result[j],string(output,"_fig_",j,".png"))
    end
    for j = 15:27
        curr_output = string(output,"_mat_",j-14,".txt")
        open(curr_output,"w") do io
            writedlm(io,result[j])
        end
    end
    curr_output = string(output,"_xaxis.txt")
    open(curr_output,"w") do io
        writedlm(io,result[28])
    end
    curr_output = string(output,"_yaxis.txt")
    open(curr_output,"w") do io
        writedlm(io,result[29])
    end
end

#takes an NxN matrix and converts NaNs to zeros
function NaNtoZero(matrix)
    N = size(matrix)[1]
    result = zeros(Float64,N,N)
    setVal = 0
    #find some value in matrix
    for i in 1:N
        for j in 1:N
            if isnan(matrix[i,j]) == false
                setVal = matrix[i,j]
                break
            end
        end
    end
    
    for i in 1:N
        for j in 1:N
            if isnan(matrix[i,j])
                result[i,j] = setVal
            else
                result[i,j] = matrix[i,j]
            end
        end
    end
    return result
end



#find information about mode projections:
#RMS of dominant mode, mode number, eigenvalue
function find_mode_info(sample)
    #project each mode onto sample
    #(except 0 mode corresponding to volume shift)
    max_rms = 0
    max_ind = 0
    largest_proj = NaN
    for i = Ne-11:size(elastic_modes)[2]-1
        time_series = zeros(length(sample))
        for j = 1:length(time_series)
            signal = sample[j][1:Ne]
            signal = signal .- mean(signal)
            norm = sqrt(sum(signal.^2))
            time_series[j] = abs(dot(signal,modes[:,i])) / (norm * sqrt(sum(modes[:,i].^2)))
        end
        rms = mean(time_series)
        if rms > max_rms
            max_rms = rms 
            max_ind = i 
        end
        if i == size(elastic_modes)[2]-1
            largest_proj = mean(time_series)
        end
    end
   
    return max_rms, max_ind, largest_proj
end