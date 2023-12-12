# Given some state u, h finds set of diff eq. of system
function h(du, u, p, t)
    # Failsafe
    if sum(abs.(u) .> 5*maximum(V0)) > 0 
        println(t)
        throw(string("Solution unstable at " ,t))
    end

    V = u[1:Ne] #volume 
    C = u[Ne+1:2Ne] #solute amount
    radius_sq = (V) ./ (pi .* Len)
    R = sqrt.(radius_sq)

    ce = C ./ V #concentration 
    ϵ = (R .- R0)./R0 #percent change from equilibrium volume

    #Volume dynamics
                  #  restoring force      nonlinear fudge    actomyosin contractions 
    r = Mq_inv \ (-(k .*ϵ) .- (kappa .*(ϵ).^3) .- ( α .*(C./C0).*(1 .- ϵ/ϵ_s)))
    μ = sum(r) / sum_Minv1
    Vdot = r .- (Minv1 .* μ)

    #Flows 
    Qmid = mid_flows(Vdot)
    Qin, Qout = in_out_flows(Qmid, Vdot)

    #Taylor dispersion
    Deff_A = Deff .* pi .* radius_sq
    DeffIN = Deff_A .*(1 .+(Qin.^2)) ./ (48 .*pi^2 .*radius_sq .*Dt.^2)
    DeffOUT = Deff_A .*(1 .+(Qout.^2)) ./ (48 .*pi^2 .*radius_sq .*Dt.^2)

    #Solve for node concentrations algebriacly
    Jin, Jout, cn, BDeff = solute_in_out_and_node_conc(Qin, Qout, ce, DeffOUT, DeffIN)

    #Concentration dynamics
            #advection                production              decay              diffusion 
    Cdot = Jin .- Jout .+ (surface_area .* β .*((1 .+ ϵ/ϵ_c) .- ce)) .+ (BDeff *cn) .- (DeffIN .+ DeffOUT).*ce

    #Save dynamics
    du .= [Vdot; Cdot]
end


#Find flow velocities given some state u
function findFlows(u)
    # Failsafe
    if  sum(abs.(u) .> 5*maximum(V0)) > 0 
        println(t)
        throw(string("Solution diverges at " ,t))
    end

    V = u[1:Ne] #volume 
    C = u[Ne+1:2Ne] #solute amount
    radius_sq = (V) ./ (pi .* Len)
    R = sqrt.(radius_sq)

    ϵ = (R .- R0)./R0 #percent change from equilibrium volume

    #Volume dynamics
                  #  restoring force      nonlinear fudge    actomyosin contractions 
    r = Mq_inv \ (-(k .*ϵ) .- (kappa .*(ϵ).^3) .- ( α .*(C./C0).*(1 .- ϵ/ϵ_s)))
    μ = sum(r) / sum_Minv1
    Vdot = r .- (Minv1 .* μ)

    #Flows 
    Qmid = mid_flows(Vdot)
    return 2*Qmid./(pi.*R.^2)
end