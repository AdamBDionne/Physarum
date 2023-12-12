using StatsBase
using LightGraphs 
using DifferentialEquations
using Plots
using Images

file_suffix = "sample_spring_seven"
p = [0.38 0.4]
fps = 3

# file_suffix = "sample_spring_six"
# p = [0.38 0.6]
#u0 = readdlm("frame_model2_u0.txt")

#sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
#time_series_rms, time_series_rms_model, avg_amp, avg_amp_model, ind_real, ind_model, indmaxm, indmaxr = project_real_and_model2(file_suffix, sol, elastic_modes, Ne, fps)
@time plot_projections2(file_suffix, time_series_rms, time_series_rms_model, avg_amp, avg_amp_model, ind_real, ind_model, indmaxm, indmaxr)