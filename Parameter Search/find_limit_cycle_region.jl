using LightGraphs 

# include("fxns/compute_jacobian.jl")
# include("fxns/model.jl")
# include("fxns/model_helper_functions.jl")

#Physarum 
    g, avg_radii, Len = import_MATLAB_network("sample_spring_nine", false)
    params = graph_precomputation(g, avg_radii, Len)

#Example graphs
    # n = 21
    # g = LightGraphs.grid([n, 1]; periodic=false)
    # params = graph_precomputation(g, [], [])

# Save parameters
Ne, Nv, avg_radii, Len, V0, kappa, α, D, β, Deff, Dt, τ, C0, E, B, K, k, L, Ldag, Mq_inv, Minv1, sum_Minv1, elastic_modes, elastic_eigvals, edgelist, out_neighbors, in_neighbors, Ein, Eout, EinT, EoutT, surface_area, radius_sq, Char_L, Char_t, R0, Cq, mu, b = params
                                    #ϵ_s     #ϵ_c     #size
#@time result = constructJacobian(g, 0.28, 0.45, 0.0, 1.5, 15)

p = [0.38 0.3] 
ϵ_s, ϵ_c = p

(~,Jv) = Zygote.forward_jacobian(u -> f(u,p[1],p[2]), [V0; V0])
            
eigs = eigvals(Jv)
num_eigs = sum((abs.(imag.(eigs))) .> 0)
imags = imag.(eigs)
real_of_imags = real.(eigs)
deleteat!(real_of_imags, abs.(imags) .< 0.001)
if num_eigs == 2*Ne-2
    if maximum(real_of_imags) > 0
        println("!")
    else
        println("BAD")
    end
else 
    println("BAD")
end