# Adam Dionne
# This script saves the data collected from a simulation
function save_sim(destination)

    odeSol = [sol[i] for i = 1:size(sol,2)]

    output3 = string(destination, "_ODE_SOL_2.txt")
    open(output3,"w") do io
        writedlm(io,odeSol)
    end

end

