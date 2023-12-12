using FFTW

# Compares model to observed time period & volumetric flows 
file_suffix = "sample_two"
source = string("Raw data/",file_suffix,"/",file_suffix)
destination = string("Processed Data/",file_suffix,"/",file_suffix)

# compare radii vs time
edge_num = 120

#model
radii_time_model = zeros(length(sol.t), Ne)
for i = 1:length(sol.t)
    radii_time_model[i,:] = sol.u[i][1:Ne] * Char_L
end

filename_radii = string(source,"_radii_time.txt")
radii_data = readdlm(filename_radii) #.* 2 * 1.83 * 10^(-6)
radii_time_real = zeros(length(split(radii_data[3],",")),Ne)

for i = 1:Ne
    broken_up = split(radii_data[i],",")
    radii_data_curr = parse.(Float64,broken_up)
    radii_time_real[:,i] = radii_data_curr * 2 * 1.83 * 10^(-6)
end

equil_radii = zeros(Ne,2)
for i = 1:Ne
    equil_radii[i,1] = i 
    equil_radii[i,2] = mean(radii_time_real[:,i])
end
equil_radii = equil_radii[sortperm(equil_radii[:,2]),:]
equil_radii_order = Int.(equil_radii[:,1])

#radii_time_model = radii_time_model[:,equil_radii_order]
#radii_time_real = radii_time_real[:,equil_radii_order]

#observed time period 
time_periods = zeros(Ne,2)
for i = 1:Ne
    real_signal = radii_time_real[:,i]
    model_signal = radii_time_model[:,i]

    #center both around 1
    real_signal = real_signal .- mean(real_signal)
    model_signal = model_signal .- mean(model_signal)

    #fourier transforms 
    F = fft(model_signal) |> fftshift
    freqs = fftfreq(length(sol.t), 10/fps) |> fftshift

    zero_ind = argmax(freqs .== 0)
    F_pos = F[zero_ind:end]
    freqs_pos = freqs[zero_ind:end]

    #time period
    maxima = partialsortperm(abs.(F_pos),2:3,rev=true)
    first_harmonic_freq = freqs_pos[maxima[1]]
    model_time_period = 1/first_harmonic_freq

    F = fft(real_signal) |> fftshift
    freqs = fftfreq(length(real_signal), 6) |> fftshift

    zero_ind = argmax(freqs .== 0)
    F_pos = F[zero_ind:end]
    freqs_pos = freqs[zero_ind:end]

    #time period
    maxima = partialsortperm(abs.(F_pos),2:3,rev=true)
    first_harmonic_freq = freqs_pos[maxima[1]]
    real_time_period = 1/first_harmonic_freq

    time_periods[i,1] = real_time_period
    time_periods[i,2] = model_time_period
end

#now find volumetric flows for real and model 

dV_real = zeros(size(radii_time_real,1),Ne)

#real
for i = 1:Ne
    #find flows 
    radii_time_curr = radii_time_real[:,i]

    #find dV
    volume_time = pi .* radii_time_curr.^2 .* Len[i]
    dV_real[2:end-1,i] = (volume_time[3:end] - volume_time[1:end-2])/(2*5)
    dV_real[1,i] = (volume_time[2] - volume_time[1]) / 5
    dV_real[end,i] = (volume_time[end] - volume_time[end-1]) / 5
end

global avg_flow_real = zeros(Ne)
for i = 1:size(radii_time_real,1)
    global avg_flow_real
    Qmid = mid_flows(dV_real[i,:]) ./ (pi .* radii_time_real[i,:].^2)
    avg_flow_real = avg_flow_real .+ Qmid.^2
end
avg_flow_real = sqrt.(avg_flow_real ./ size(radii_time_real,1))

#model
flows = [findFlows(sol(t_i+(i/fps))) ./ (sol(t_i+(i/fps))[1:Ne] ./ Len) for i in 1:fps*(t_f-t_i)]
global avg_flow_model = zeros(Ne)
for i = 1:length(flows)
    global avg_flow_model
    avg_flow_model += flows[i][equil_radii_order].^2
end
avg_flow_model = sqrt.(avg_flow_model ./ length(flows))



output = string("vel_model.txt")
open(output,"w") do io
    writedlm(io,avg_flow_model)
end

output = string("vel_real.txt")
open(output,"w") do io
    writedlm(io,avg_flow_real)
end



