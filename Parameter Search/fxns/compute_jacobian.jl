using FFTW
using DifferentialEquations
using LinearAlgebra
using Statistics
using LightGraphs
using Zygote
using Plots
using Images
using Combinatorics

include("hilbert.jl")
include("phase_diagram_helper_functions.jl")
#include("model.jl")

#scan over parameter space, and see what the Jacobian is telling us 
function constructJacobian(g,f1Init,f1Final,f2Init,f2Final,size_grid)

    #initialize matrices that store diagram data 
    phaseMat1 = zeros(Float64,size_grid,size_grid)
    phaseMat2 = zeros(Float64,size_grid,size_grid)

    #step values in parameter space
    f1Step = (f1Final - f1Init) / size_grid
    f2Step = (f2Final - f2Init) / size_grid

    #iterate through each pixel wanted
    for i in 1:size_grid
        #update progress
        print(string(i," "))
        for j in 1:size_grid
            #determine parameters 
            p = [(f1Init + f1Step*i) (f2Init + f2Step*j)]
            ϵ_s, ϵ_c = p
            #compute Jacobian
            (~,Jv) = Zygote.forward_jacobian(u -> f(u,p[1],p[2]), [V0; C0])
            
            eigs = eigvals(Jv)
            num_eigs = sum((abs.(imag.(eigs))) .> 0)
            phaseMat2[i,j] = num_eigs
            imags = imag.(eigs)
            real_of_imags = real.(eigs)
            deleteat!(real_of_imags, abs.(imags) .< 0.001)
            if isempty(real_of_imags)
                phaseMat2[i,j] = NaN
            else
                phaseMat2[i,j] = maximum(real_of_imags)
            end
            if num_eigs == 2*Ne-2
                if maximum(real_of_imags) > 0
                    phaseMat1[i,j] = 1
                end
                # if maximum(real.(eigs.*(abs.(imag.(eigs)) .> 0))) > 0
                #     phaseMat1[i,j] = 1
                # end
            else
                phaseMat1[i,j] = NaN 
            end
        end
    end
    xaxis = [f2Init+f2Step*i for i in 1:size_grid]
    yaxis = [f1Init+f1Step*i for i in 1:size_grid]
    #heatmap plot

    #scales
    colorscale = cgrad(:roma, scale = :log)


    #find clim
    zerosMat1 = NaNtoZero(phaseMat1)
    zerosMat2 = NaNtoZero(phaseMat2)


    p1 = heatmap(xaxis,yaxis,phaseMat1,title="Phase Diagram Comparing ϵs and ϵc.",xlabel="ϵc",ylabel="ϵs",cbar=true, seriescolor = colorscale, clim = (minimum(zerosMat1),maximum(zerosMat1)), colorbar_title = "Limit Cycle (Yay/Nay)")
    p2 = heatmap(xaxis,yaxis,phaseMat2,title="Phase Diagram Comparing ϵs and ϵc.",xlabel="ϵc",ylabel="ϵs",cbar=true, seriescolor = colorscale, clim = (minimum(zerosMat2),maximum(zerosMat2)), colorbar_title = "# Complex Eigs")
    return p1, p2, phaseMat1, phaseMat2
end