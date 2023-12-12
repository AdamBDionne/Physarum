using LightGraphs
using DifferentialEquations
using Statistics
using Images
using StatsBase

# file_suffix = "sample_spring_one"
# p = [0.38 0.8]
# fps = 5

# sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
# @time time_series_rms, time_series_rms_model, avg_amp, avg_amp_model = project_real_and_model(file_suffix, sol, elastic_modes, Ne, fps)
# @time plot_projections(file_suffix, time_series_rms, time_series_rms_model, avg_amp, avg_amp_model)



file_suffix = "sample_spring_two"
p = [0.38 0.8]
fps = 3

sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
@time time_series_rms, time_series_rms_model, avg_amp, avg_amp_model = project_real_and_model(file_suffix, sol, elastic_modes, Ne, fps)
@time plot_projections(file_suffix, time_series_rms, time_series_rms_model, avg_amp, avg_amp_model)

file_suffix = "sample_spring_three"
p = [0.38 0.4]
fps = 3

sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
@time time_series_rms, time_series_rms_model, avg_amp, avg_amp_model = project_real_and_model(file_suffix, sol, elastic_modes, Ne, fps)
@time plot_projections(file_suffix, time_series_rms, time_series_rms_model, avg_amp, avg_amp_model)

file_suffix = "sample_spring_four"
p = [0.38 0.3]
fps = 3

sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
@time time_series_rms, time_series_rms_model, avg_amp, avg_amp_model = project_real_and_model(file_suffix, sol, elastic_modes, Ne, fps)
@time plot_projections(file_suffix, time_series_rms, time_series_rms_model, avg_amp, avg_amp_model)

file_suffix = "sample_spring_five"
p = [0.38 0.7]
fps = 3

sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
@time time_series_rms, time_series_rms_model, avg_amp, avg_amp_model = project_real_and_model(file_suffix, sol, elastic_modes, Ne, fps)
@time plot_projections(file_suffix, time_series_rms, time_series_rms_model, avg_amp, avg_amp_model)

file_suffix = "sample_spring_six"
p = [0.38 0.6]
fps = 3

sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
@time time_series_rms, time_series_rms_model, avg_amp, avg_amp_model = project_real_and_model(file_suffix, sol, elastic_modes, Ne, fps)
@time plot_projections(file_suffix, time_series_rms, time_series_rms_model, avg_amp, avg_amp_model)

file_suffix = "sample_spring_seven"
p = [0.38 0.4]
fps = 3

sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
@time time_series_rms, time_series_rms_model, avg_amp, avg_amp_model = project_real_and_model(file_suffix, sol, elastic_modes, Ne, fps)
@time plot_projections(file_suffix, time_series_rms, time_series_rms_model, avg_amp, avg_amp_model)

file_suffix = "sample_spring_eight"
p = [0.38 0.2]
fps = 3

sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
@time time_series_rms, time_series_rms_model, avg_amp, avg_amp_model = project_real_and_model(file_suffix, sol, elastic_modes, Ne, fps)
@time plot_projections(file_suffix, time_series_rms, time_series_rms_model, avg_amp, avg_amp_model)

file_suffix = "sample_spring_nine"
p = [0.38 0.3]
fps = 3

sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
@time time_series_rms, time_series_rms_model, avg_amp, avg_amp_model = project_real_and_model(file_suffix, sol, elastic_modes, Ne, fps)
@time plot_projections(file_suffix, time_series_rms, time_series_rms_model, avg_amp, avg_amp_model)


