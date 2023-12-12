using LinearAlgebra
using Statistics
using DelimitedFiles

# Defines all neccessary parameters for model
function graph_precomputation(g, avg_radii, Len)
    #Store number of edges and vertices
    Ne, Nv = ne(g), nv(g)

    # Check if radii or length are pre-set. If not, give default values 
    if isempty(avg_radii)
        avg_radii = ones(Ne) .* 46 * 1.83 * 10^(-6);
    end

    if isempty(Len)
        Len = ones(Ne) .* 46 * 5 .* 1.83 * 10^(-6);
    end

    #Parameter Values (unscaled)
    V0 = pi .* (avg_radii).^2 .* Len
    mu = 1.5 * 10^(-3) # Pa*s viscosity   SET
    R = 46 * 1.83 * 10^(-6) #meters     
    #Lp = 46* 5 * 1.83 * 10^(-6) #meters  
    kp = 10 #pascals    SET
    eta = 240 #pascals seconds    SET
    #b = tau / R #damping costnat kg s^(-1)
    kappa = 2000.0 # (Pascals) (kappa) nonlinear elastic constant
    α = 30.0 # (Pascals) active concentration dynamics parameter
    D = 3.33 * 10^(-10) #m^2 / s
    β = avg_radii ./ (2*96) # (L/s)  production scale 

    # Set characteristic scales 
    #Char_V = (pi*R^2)*Lp # m^3
    Char_V = mean(V0)
    Char_L = (Char_V)^(1/3) # m 
    Char_t = 15 # s 
    # note: can set pressure scale independantly through mass scaling
    Char_p = 100 # Pa 



    # Nondimensionalize parameters
    mu = mu * (1/(Char_p*Char_t)) # Pa*s viscosity 
    R = R * (1/Char_L) #meters
    k = kp * (1/Char_p) #pascals
    b = eta * (1/(Char_p*Char_t)) #pascals*s
   # b = tau / R
    kappa = kappa * (1/Char_p) # (Pascals) (kappa) nonlinear elastic constant
    α = α * (1/Char_p) # (Pascals) active concentration dynamics parameter
    Dt = D * (Char_t/(Char_L^2)) # m^2 / s
    β = β * (Char_t/Char_L) # (1/s)  production scale 

    avg_radii = avg_radii .* (1/Char_L)
    #b = tau ./ avg_radii
    Len = Len .* (1/Char_L)
    V0 = vec(V0 .* (1/Char_V))
    C0 = V0
    Deff = (2 .* Dt) ./ (Len)

        K = (pi.*avg_radii.^4)./(8*mu.*Len)
        τ = diagm(0 => vec((b ./(4 .* pi .* avg_radii.^2 .* Len)) .+ (1 ./ (12*K))))


    radius_sq = V0 ./ (pi.*Len)
    surface_area = 2 .*pi .*sqrt.(radius_sq) .* Len

    ## Pre-computation for model.

    #Here we calculate the Laplacian L and 
    #matrix M - see Henrik's notes eqn (51)
    E = incidence_matrix(g; oriented=true)
    B = incidence_matrix(g; oriented=false)

    L = -E*(diagm(0 => vec(K)))*E'
    Ldag = pinv(Array(L))
    
    Cq = -0.25B'*Ldag*B
    Mq = τ - 0.25B'*Ldag*B

        # Pre-calculation for simulation speed ups
        Mq_inv = factorize(Mq)
        Minv1 = Mq_inv \ ones(Ne)
        sum_Minv1 = sum(Minv1)

    # Modes
    dim = size(Mq)[1]
    scaling = diagm(0 => vec(2 .* pi .* Len .* avg_radii.^2))
    A = -k*(I-((ones(dim)*ones(dim)'*inv(Mq)/(ones(dim)'*inv(Mq)*ones(dim)))))
    elastic_M = eigen(A,Mq*scaling)
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

    R0 = sqrt.(V0 ./ (pi .* Len))

    return Ne, Nv, avg_radii, Len, V0, kappa, α, D, β, Deff, Dt, τ, C0, E, B, K, k, L, Ldag, Mq_inv, Minv1, sum_Minv1, elastic_modes, elastic_eigvals, edgelist, out_neighbors, in_neighbors, Ein, Eout, EinT, EoutT, surface_area, radius_sq, Char_L, Char_t, R0, Cq, mu, b
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


#returns graph, lengths, radii, 
#and can save in standard format 
function import_MATLAB_network(file_suffix, save_result)
    source = string("Data/",file_suffix,"/",file_suffix)
    #destination = string("Processed data/",file_suffix,"/",file_suffix)
    #helper function to find what node index corresponds to what edge
    function find_node_ind(x,y,node_data)
        for i = 1:size(node_data)[1]
            if node_data[i,1] == x && node_data[i,2] == y
                return i
            end
        end
        return 0
    end

    #data to be read in
    filename_edges = string(source,"_network_edges.csv")
    filename_nodes = string(source,"_network_nodes.txt")

    edge_data = readdlm(filename_edges)
    node_data = readdlm(filename_nodes)
    node_data = round.(Int,node_data)


    #build graph
    g = SimpleGraph(size(node_data)[1])

    pix_len_order = zeros(size(edge_data)[1]-1)
    pix_len_order2 = zeros(size(edge_data)[1]-1)
  
    edge_nodes = []
    #add in edges
    for i = 2:size(edge_data)[1]
        #get edge coords
        curr_edge = edge_data[i];
        curr_edge_dat = split(curr_edge,',')
        edge_pieces = parse.(Int,curr_edge_dat[1:4])
        #find corresponding node indices
        #node_one = [edge_pieces[1] edge_pieces[2]]
        #node_two = [edge_pieces[3] edge_pieces[4]]
        node_one = find_node_ind(edge_pieces[1],edge_pieces[2],node_data)
        node_two = find_node_ind(edge_pieces[3],edge_pieces[4],node_data)
        push!(edge_nodes, [node_one; node_two])
        #edge_nodes[i-1,1] = node_one 
        #edge_nodes[i-1,2] = node_two
    
        add_edge!(g,node_one,node_two)
        pix_len_order[i-1] = parse(Float64,curr_edge_dat[5])
        pix_len_order2[i-1] = ceil(sqrt((edge_pieces[1] - edge_pieces[3])^2 + (edge_pieces[2] - edge_pieces[4])^2))
    end

    #also load in average radii 
    filename_radii = string(source,"_radii_time.txt")
    radii_data = readdlm(filename_radii) #.* 1.83 * 10^(-6)
    avg_radii = zeros(ne(g))

    for j = 1:ne(g)
        broken_up = split(radii_data[j],",")
        radii_data_curr = parse.(Float64,broken_up)
        avg_radii[j] = mean(radii_data_curr) * 1.83 * 10^(-6)

    end
     #avg_radii = []
    #get lengths 
    Len = ones(ne(g))
    for (i,e) in enumerate(edges(g))
        Len[i] = pix_len_order[i] * 1.83*10^(-6)
    end

    if save_result
        #save a translation for graph structures from julia to matlab
        #(the ordering is different)
        translator = zeros(ne(g))
        for (i,e) in enumerate(edges(g))
            t = findall(x->x==[src(e),dst(e)],edge_nodes)
            if t == []
                translator[i] = -findall(x->x==[dst(e),src(e)],edge_nodes)[1]
            else
                translator[i] = t[1]
            end
        end
        output3 = string(source,"_translator.txt")
        open(output3,"w") do io
            writedlm(io,translator)
        end

        #now save graph in the same format as used for our animations..
        edgesM = zeros(Int,ne(g),3)
        for (i,e) in enumerate(edges(g))
            edgesM[i,1] = src(e)
            edgesM[i,2] = dst(e)
            edgesM[i,3] = pix_len_order2[Int(abs(translator[i]))]
        end

        #save
        output1 = string(source,"_EDGES.txt")
        open(output1,"w") do io
            writedlm(io,edgesM)
        end

        #similarly save nodes...
        nodesM = zeros(Int,size(node_data)[1],3)
        for i = 1:size(node_data)[1]
            nodesM[i,1] = i
            nodesM[i,2] = node_data[i,1]
            nodesM[i,3] = node_data[i,2]
        end

        #save
        output2 = string(source, "_NODES.txt")
        open(output2,"w") do io
            writedlm(io,nodesM)
        end

        
    end

    #create a translation key that simply rearranges
    #rearrange the graph structure in matlab? 

    return g, avg_radii, Len
end




    
