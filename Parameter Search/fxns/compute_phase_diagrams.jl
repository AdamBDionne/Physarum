#--- packages
using FFTW
using DifferentialEquations
using LinearAlgebra
using Statistics
using LightGraphs
using Zygote
using Plots
using Images
using Combinatorics
using StatsBase

include("hilbert.jl")
include("phase_diagram_helper_functions.jl")
include("model.jl")

function constructPhaseDiagram(g,f1Init,f1Final,f2Init,f2Final,size_grid)

    #initialize matrices that store diagram data 
    phaseMat1 = zeros(Float64,size_grid,size_grid)   #time period 
    phaseMat2 = zeros(Float64,size_grid,size_grid)   #amplitude
    phaseMat3 = zeros(Float64,size_grid,size_grid)   #flow velocity
    phaseMat4 = zeros(Float64,size_grid,size_grid)   #global order 
    phaseMat5 = zeros(Float64,size_grid,size_grid)   #local order
    phaseMat6 = zeros(Float64,size_grid,size_grid)   #largest Re of complex
    phaseMat7 = zeros(Float64,size_grid,size_grid)   #largest proj
    phaseMat8 = zeros(Float64,size_grid,size_grid)   #largest proj ind
    phaseMat9 = zeros(Float64,size_grid,size_grid)   #largest mode proj
    phaseMat10 = zeros(Float64,size_grid,size_grid)  #THDF

    #step values in parameter space
    f1Step = (f1Final - f1Init) / size_grid
    f2Step = (f2Final - f2Init) / size_grid

    #iterate through each pixel wanted
    for i in 1:size_grid
        #update progress
        print(string(i," "))
        for j in 1:size_grid
           
            #determine parameters 
            p = [(f1Init + f1Step*i); (f2Init + f2Step*j)]

            #compute Jacobian
            (~,Jv) = Zygote.forward_jacobian(u -> f(u,p[1],p[2]), [V0; V0])
                
            eigs = eigvals(Jv)
            num_eigs = sum((abs.(imag.(eigs))) .> 0)

            imags = imag.(eigs)
            real_of_imags = real.(eigs)
            deleteat!(real_of_imags, abs.(imags) .< 0.001)

            limit_cycle = false
            if isempty(real_of_imags)
                phaseMat6[i,j] = NaN
            else
                max_real = maximum(real_of_imags)
                phaseMat6[i,j] = max_real
                if max_real > 0 && (num_eigs == 2Ne-2)
                    limit_cycle = true
                end
            end

            #check if in limit cycle...
            #(there are always two zero eigenvalues for our jacobians)
            if limit_cycle
                #solve the system
                #try
                    T = 200
                    tspan = (0.0,T)
                    shift = elastic_modes[:,Ne-1] - elastic_modes[:,Ne-2]
                    u0 = [V0; C0 .+ (0.05 .*(shift ./ maximum(abs.(shift))))]
                    ϵ_s, ϵ_c = p
                    prob = ODEProblem(h,u0,tspan,p)
                    sol = solve(prob, RK4(),saveat=50:0.01:200)
                    sampleTime = 0.01
                    T0 = 50
                    Tf = 200
                    sample = sol(T0:sampleTime:Tf).u
                    
                    # t_period, amp, THDF = fft_info(sample, sampleTime)
                    # phaseMat1[i,j] = t_period
                    # phaseMat2[i,j] = amp
                    # phaseMat10[i,j] = THDF

                    # global_order, local_order = find_order(sample)
                    # phaseMat4[i,j] = global_order
                    # phaseMat5[i,j] = local_order

                    # flow_amp = find_flow_amp(sample, sampleTime, p)
                    # phaseMat3[i,j] = flow_amp

                    # max_rms, max_ind, largest_proj = find_mode_info(sample)
                    # phaseMat7[i,j] = max_rms
                    # phaseMat8[i,j] = max_ind
                    # phaseMat9[i,j] = largest_proj

                    #set amount simulated
                    # global num_parts = 3*Ne
                    # global num_per = 3
                    # global dt = 2/sampleTime

                    # #set initial conditions
                    # pos_0 = []
                    # for i = 1:num_per
                    #     append!(pos_0,i.*Len/(num_per+1))
                    # end

                    # ind_0 = []
                    # for i = 1:num_per
                    #     append!(ind_0,[i for i = 1:Ne])
                    # end

                    #find neighbors for every vertex
                    # global every_neighbor = []
                    # for i=1:Nv
                    #     # pos neighboring edges 
                    #     push!(every_neighbor, [in_neighbors[i] ; out_neighbors[i]])
                    # end

                    #get needed values

                    #model 
                    global Qin_t = findFlowsOut.(sol.u)
                    global Qout_t = findFlowsIn.(sol.u)
                    global A_t = []
                    for i = 1:length(sol.t)
                        push!(A_t,sol.u[i][1:Ne] ./ (Len))
                    end

                    pos_time, ind_time = simulate(A_t, Qin_t, Qout_t, pos_0, ind_0)
                    t = find_average_displacement(pos_time, ind_time)
                    phaseMat1[i,j] = mean(t)

                # catch
                #     phaseMat1[i,j] = NaN
                #     phaseMat2[i,j] = NaN
                #     phaseMat3[i,j] = NaN
                #     phaseMat4[i,j] = NaN
                #     phaseMat5[i,j] = NaN
                #     phaseMat7[i,j] = NaN
                #     phaseMat8[i,j] = NaN
                #     phaseMat9[i,j] = NaN
                #     phaseMat10[i,j] = NaN
                # end
            else
                phaseMat1[i,j] = NaN
                phaseMat2[i,j] = NaN
                phaseMat3[i,j] = NaN
                phaseMat4[i,j] = NaN
                phaseMat5[i,j] = NaN
                phaseMat7[i,j] = NaN
                phaseMat8[i,j] = NaN
                phaseMat9[i,j] = NaN
                phaseMat10[i,j] = NaN
            end
        end
    end
    xaxis = [f2Init+f2Step*i for i in 1:size_grid]
    yaxis = [f1Init+f1Step*i for i in 1:size_grid]
    #heatmap plot

    #scales
        colorscale = cgrad(:roma, scale = :log)
        colorscale2 = cgrad(:thermal)
        colorscale3 = cgrad(:cyclic_wrwbw_40_90_c42_n256)

    #find clim
    zerosMat1 = NaNtoZero(phaseMat1)
    zerosMat2 = NaNtoZero(phaseMat2)
    zerosMat3 = NaNtoZero(phaseMat3)
    zerosMat4 = NaNtoZero(phaseMat4)
    zerosMat5 = NaNtoZero(phaseMat5)
    zerosMat6 = NaNtoZero(phaseMat6)
    zerosMat7 = NaNtoZero(phaseMat7)
    zerosMat8 = NaNtoZero(phaseMat8)
    zerosMat9 = NaNtoZero(phaseMat9)
    zerosMat10 = NaNtoZero(phaseMat10)

    p1 = heatmap(xaxis,yaxis,phaseMat1,title="Phase Diagram Comparing ϵs and ϵc.",xlabel="ϵc",ylabel="ϵs",cbar=true, seriescolor = colorscale2, clim = (minimum(zerosMat1),maximum(zerosMat1)), colorbar_title = "Time Period")
    p2 = heatmap(xaxis,yaxis,phaseMat2,title="Phase Diagram Comparing ϵs and ϵc.",xlabel="ϵc",ylabel="ϵs",cbar=true, seriescolor = colorscale2, clim = (minimum(zerosMat2),maximum(zerosMat2)), colorbar_title = "Amplitude")
    p3 = heatmap(xaxis,yaxis,phaseMat3,title="Phase Diagram Comparing ϵs and ϵc.",xlabel="ϵc",ylabel="ϵs",cbar=true, seriescolor = colorscale2, clim = (minimum(zerosMat3),maximum(zerosMat3)), colorbar_title = "Flow Velocity Amplitude")
    p4 = heatmap(xaxis,yaxis,phaseMat4,title="Phase Diagram Comparing ϵs and ϵc.",xlabel="ϵc",ylabel="ϵs",cbar=true, seriescolor = colorscale2, clim = (minimum(zerosMat4),maximum(zerosMat4)), colorbar_title = "Global Order")
    p5 = heatmap(xaxis,yaxis,phaseMat5,title="Phase Diagram Comparing ϵs and ϵc.",xlabel="ϵc",ylabel="ϵs",cbar=true, seriescolor = colorscale2, clim = (minimum(zerosMat5),maximum(zerosMat5)), colorbar_title = "Local Order")
    p6 = heatmap(xaxis,yaxis,phaseMat6,title="Phase Diagram Comparing ϵs and ϵc.",xlabel="ϵc",ylabel="ϵs",cbar=true, seriescolor = colorscale2, clim = (minimum(zerosMat6),maximum(zerosMat6)), colorbar_title = "Largest Real Part Among Complex Eigs")
    p7 = heatmap(xaxis,yaxis,phaseMat7,title="Phase Diagram Comparing ϵs and ϵc.",xlabel="ϵc",ylabel="ϵs",cbar=true, seriescolor = colorscale2, clim = (minimum(zerosMat7),maximum(zerosMat7)), colorbar_title = "Largest Projection RMS")
    p8 = heatmap(xaxis,yaxis,phaseMat8,title="Phase Diagram Comparing ϵs and ϵc.",xlabel="ϵc",ylabel="ϵs",cbar=true, seriescolor = colorscale2, clim = (minimum(zerosMat8),maximum(zerosMat8)), colorbar_title = "Largest Projection Index")
    p9 = heatmap(xaxis,yaxis,phaseMat9,title="Phase Diagram Comparing ϵs and ϵc.",xlabel="ϵc",ylabel="ϵs",cbar=true, seriescolor = colorscale2, clim = (minimum(zerosMat9),maximum(zerosMat9)), colorbar_title = "Largest Excitation Projection")
    p10 = heatmap(xaxis,yaxis,phaseMat10,title="Phase Diagram Comparing ϵs and ϵc.",xlabel="ϵc",ylabel="ϵs",cbar=true, seriescolor = colorscale2, clim = (minimum(zerosMat10),maximum(zerosMat10)), colorbar_title = "THDF")

    return p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, phaseMat1, phaseMat2, phaseMat3, phaseMat4, phaseMat5, phaseMat6, phaseMat7, phaseMat8, phaseMat9, phaseMat10, xaxis, yaxis
end