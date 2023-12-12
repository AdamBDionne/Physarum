#function that defines set of diff eq. for system 
#given initial condition u 
function h(du, u, p, t)
    ϵ_s, ϵ_c = p
    # Failsafe: if we are in the wrong regime, the solution 
    # will diverge. This might cause the program to halt,  
    # so we throw an error if any vessel volume is too large
    if  sum(abs.(u) .> 20*maximum(V0)) > 0 
        println(t)
        throw(string("Solution diverges at " ,t))
    end

    V = u[1:Ne] #volume 
    C = u[Ne+1:2Ne] #solute amount
    radius_sq = (V) ./ (pi .* Len)
    R = sqrt.(radius_sq)

    ce = C ./ V #concentration 
    ϵ = (R .- R0)./R0 #percent change from equilibrium volume

    
    # nonlinear conductivities and damping:
        # τ = diagm(0 => vec((b ./(4 .* pi^2 .* radius_sq .* Len.^2))))
        # Mq = τ + Cq
        # Mq_inv = factorize(Mq)
        # Minv1 = Mq_inv \ ones(Ne)
        # sum_Minv1 = sum(Minv1)

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
            #advection                production   decay              diffusion 
    Cdot = Jin .- Jout .+ (surface_area .* β .*((1 .+ ϵ/ϵ_c) .- ce)) .+ (BDeff *cn) .- (DeffIN .+ DeffOUT).*ce

    #Save dynamics
    du .= [Vdot; Cdot]
end


#function that defines set of diff eq. for system 
#given initial condition u 
function f(u, ϵ_s, ϵ_c)
         # Failsafe: if we are in the wrong regime, the solution 
    # will diverge. This might cause the program to halt,  
    # so we throw an error if any vessel volume is too large
    if  sum(abs.(u) .> 5*maximum(V0)) > 0 
        println(t)
        throw(string("Solution diverges at " ,t))
    end

    V = u[1:Ne] #volume 
    C = u[Ne+1:2Ne] #solute amount
    radius_sq = (V) ./ (pi .* Len)
    R = sqrt.(radius_sq)

    ce = C ./ V #concentration 
    ϵ = (R .- R0)./R0 #percent change from equilibrium volume

    
    # nonlinear conductivities and damping:
        #K = (pi.*R.^4)./(8*mu.*Len)
        # τ = diagm(0 => vec((b ./(4 .* pi .* R .* Len .* R0)) ))#.+ (1 ./ (12*K))))
        # Mq = τ + Cq
        # Mq_inv = factorize(Mq)
        # Minv1 = Mq_inv \ ones(Ne)
        # sum_Minv1 = sum(Minv1)

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
            #advection                production   decay              diffusion 
    Cdot = Jin .- Jout .+ (surface_area .* β .*((1 .+ ϵ/ϵ_c) .- ce)) .+ (BDeff *cn) .- (DeffIN .+ DeffOUT).*ce

    return [Vdot; Cdot]
end