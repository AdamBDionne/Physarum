#function that defines set of diff eq. for system 
#given initial condition u 
function h(du, u, p, t)
    # Failsafe: if we are in the wrong regime, the solution 
    # will diverge. This might cause the program to halt, 
    # so we throw an error if any vessel volume is too large
    if  sum(abs.(u) .> 6) > 0 
        println(t)
        throw(string("Solution diverges at " ,t))
    end
  

    V = u[1:Ne] #volume 
    Vdot = du[1:Ne] #change in volume 
    C = u[Ne+1:2Ne] #solute amount 
    Cdot = du[Ne+1:2Ne] #change in solute amount 

    ce = C ./ V #concentration 
    ϵ = (V .- V0)./V0 #percent change from equilibrium volume

    #Volume dynamics
                  #  restoring force      nonlinear fudge                 actomyosin contractions 
    r = Mq_inv \ (-(k .*ϵ .*V0) .- (kappa .* V0 .*(ϵ).^3) .- (surface_area .* α .*(C./C0).*(1 .- ϵ/ϵ_s)))
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
            #advection       production (stretch activated)   decay              diffusion 
    Cdot = Jin .- Jout .+ (surface_area .* β .*((1 .+ ϵ/ϵ_c) .- ce)) .+ (BDeff *cn) .- (DeffIN .+ DeffOUT).*ce

    #Save dynamics
    du .= [Vdot; Cdot]
end

