using LinearAlgebra
using Statistics

# Defines all neccessary parameters for model
function graph_precomputation(g, avg_radii, Len)
    #Store number of edges and vertices
    Ne, Nv = ne(g), nv(g)

    # Check if radii or length are pre-set. If not, give default values 
    if isempty(avg_radii)
        avg_radii = ones(Ne) .* 30 * 1.83 * 10^(-6);
    end

    if isempty(Len)
        Len = ones(Ne) .* 120 .* 1.83 * 10^(-6);
    end

    #Parameter Values (unscaled)
    V0 = pi .* (avg_radii).^2 .* Len
    mu = 1.5 * 10^(-3) # Pa*s viscosity
    R = 30 * 1.83 * 10^(-6) #meters
    Lp = 120 * 1.83 * 10^(-6) #meters  
    kp = 10 #pascals
    tau = 240 #pascals seconds
    kappa = 1000.0 # (Pascals) (kappa) nonlinear elastic constant
    α = 30.0 # (Pascals) active concentration dynamics parameter
    D = 3.33 * 10^(-10) #m^2 / s
    β = R/(2*96) # (L/s)  production scale 

    # Set characteristic scales 
    Char_V = (pi*R^2)*Lp # m^3
    Char_L = (Char_V)^(1/3) # m 
    Char_t = 10 # s 
    # note: can set pressure scale independantly through mass scaling
    Char_p = 100 # Pa 


    # Nondimensionalize parameters
    mu = mu * (1/(Char_p*Char_t)) # Pa*s viscosity
    R = R * (1/Char_L) #meters
    k = kp * (1/Char_p) #pascals
    tau = tau * (1/(Char_p*Char_t)) #pascals*s
    kappa = kappa * (1/Char_p) # (Pascals) (kappa) nonlinear elastic constant
    α = α * (1/Char_p) # (Pascals) active concentration dynamics parameter
    Dt = D * (Char_t/(Char_L^2)) # m^2 / s
    β = β * (Char_t/Char_L) # (1/s)  production scale 

    avg_radii = avg_radii .* (1/Char_L)
    Len = Len .* (1/Char_L)
    V0 = vec(V0 .* (1/Char_V))
    C0 = V0
    Deff = (2 .* Dt) ./ (Len)

        Kp_avg = (pi * mean(avg_radii)^4) / (8 * mu * mean(Len))
        b = tau - (1/(12*Kp_avg))

        K = (pi.*avg_radii.^4)./(8*mu.*Len)
        τ = diagm(0 => vec(b .+ (1 ./ (12*K))))

    ϵ_s = 0.5 
    ϵ_c = 0.06

    radius_sq = V0 ./ (pi.*Len)
    surface_area = 2 .*pi .*sqrt.(radius_sq) .* Len

    ## Pre-computation for model.

    #Here we calculate the Laplacian L and 
    #matrix M - see Henrik's notes eqn (51)
    E = incidence_matrix(g; oriented=true)
    B = incidence_matrix(g; oriented=false)

    L = -E*(diagm(0 => vec(K)))*E'
    Ldag = pinv(Array(L))

        Mq = τ - 0.25B'*Ldag*B

        # Pre-calculation for simulation speed ups
        Mq_inv = factorize(Mq)
        Minv1 = Mq_inv \ ones(Ne)
        sum_Minv1 = sum(Minv1)

    # Modes
    dim = size(Mq)[1]
    A = -k*(I-((ones(dim)*ones(dim)'*inv(Mq)/(ones(dim)'*inv(Mq)*ones(dim)))))
    elastic_M = eigen(A,Mq)
    elastic_modes = elastic_M.vectors
    elastic_eigvals = elastic_M.values

    edgelist = vcat([[src(e) dst(e)] for e in edges(g)]...)
    out_neighbors = []
    in_neighbors = []

    for i=1:Nv
        # pos neighboring edges 
        push!(out_neighbors, findall(edgelist[:,1] .== i))
        push!(in_neighbors, findall(edgelist[:,2] .== i))
    end

    # For helper functionsm logical E matrices
    Ein = E .* (E .< 0 ) .* (-1)
    Eout = E .* (E .> 0)
    EinT = E' .* (E' .< 0 ) .* (-1)
    EoutT = E' .* (E' .> 0)

    return Ne, Nv, avg_radii, Len, V0, kappa, α, D, β, Deff, Dt, τ, C0, E, B, K, k, L, Ldag, Mq_inv, Minv1, sum_Minv1, elastic_modes, elastic_eigvals, edgelist, out_neighbors, in_neighbors, Ein, Eout, EinT, EoutT, ϵ_s, ϵ_c, surface_area, radius_sq
end


# Calculate steady-state flows
function mid_flows(Vdot)
    # calculate the mid-point (average) flows from Vdot 
    Qmid = (0.5K).*(E'*Ldag*B*Vdot)
    return Qmid
end

# Calculate in and out flows
function in_out_flows(Qmid, Vdot)
    Qin = Qmid .+ 0.5Vdot
    Qout = Qmid .- 0.5Vdot
    return Qin, Qout
end

#find nodal concentrations along with
#solute in and out 
function solute_in_out_and_node_conc(Qin,Qout,ce,DeffOUT,DeffIN)
    #Need to seperate vectors and matrices based on sign 
    #to handle nuances of where concentration is flowing from 
    #(either the edge or nod depending on sign)

    #See Henrik's notes eqn (78)
    QinPOS = Qin .* (Qin .> 0)
    QinNEG = Qin .* (Qin .< 0)
    QoutPOS = Qout .* (Qout .> 0)
    QoutNEG = Qout .* (Qout .< 0)
    # Ein = E .* (E .< 0 ) .* (-1)
    # Eout = E .* (E .> 0)
    # EinT = E' .* (E' .< 0 ) .* (-1)
    # EoutT = E' .* (E' .> 0)

    #solve algebriac equation for nodal concentrations
    cNum = (Ein * ((DeffOUT .+ QoutPOS) .* ce)) .+ (Eout * ((DeffIN .- QinNEG) .* ce))
    cDen = (Ein * (DeffOUT .- QoutNEG) ) .+ (Eout * (DeffIN .+ QinPOS) ) 

    #decompose by going going into node v.s. going out of node 

    cn = cNum ./ cDen

    #find solute in and out flows
    Jin = (QinPOS .* (EoutT*cn)) .+ (QinNEG .* ce)
    Jout = (QoutPOS .* ce) .+ (QoutNEG .* (EinT*cn))

    #construct BDeff
    BDeff = (EoutT .* DeffIN) .+ (EinT .* DeffOUT)

    return Jin, Jout, cn, BDeff
end

