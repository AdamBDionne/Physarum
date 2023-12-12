using LightGraphs
using DifferentialEquations
using Statistics
using Images
using StatsBase

fps = 20
file_suffix = "sample_spring_one"
p = [0.38 0.8]


sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
save_sim(string("Raw Data/",file_suffix,"/",file_suffix))


file_suffix = "sample_spring_two"
p = [0.38 0.8]


sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
save_sim(string("Raw Data/",file_suffix,"/",file_suffix))

file_suffix = "sample_spring_three"
p = [0.38 0.4]


sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
save_sim(string("Raw Data/",file_suffix,"/",file_suffix))

file_suffix = "sample_spring_four"
p = [0.38 0.3]


sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
save_sim(string("Raw Data/",file_suffix,"/",file_suffix))

file_suffix = "sample_spring_five"
p = [0.38 0.7]


sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
save_sim(string("Raw Data/",file_suffix,"/",file_suffix))


file_suffix = "sample_spring_six"
p = [0.38 0.6]


sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
save_sim(string("Raw Data/",file_suffix,"/",file_suffix))


file_suffix = "sample_spring_seven"
p = [0.38 0.4]


sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
save_sim(string("Raw Data/",file_suffix,"/",file_suffix))


file_suffix = "sample_spring_eight"
p = [0.38 0.2]


sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
save_sim(string("Raw Data/",file_suffix,"/",file_suffix))


file_suffix = "sample_spring_nine"
p = [0.38 0.3]


sol, elastic_modes, Ne = @time run_model(file_suffix, p, fps)
save_sim(string("Raw Data/",file_suffix,"/",file_suffix))


