#########
#clean hanging nodes (visual processing artifact), and convert data structures from matlab -> julia
#########
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
filename_radii = string(source,"_radii_time.txt")

edge_data = readdlm(filename_edges)
node_data = readdlm(filename_nodes)
node_data = round.(Int,node_data)
radii_data = readdlm(filename_radii)

# Define a struct for edge information
struct EdgeInfo
    edge_idx::Int
    in_node_idx::Int
    out_node_idx::Int
    edge_length::Float64
    in_node_pos::Tuple{Int, Int}
    out_node_pos::Tuple{Int, Int}
    avg_radii::Float64
end

# Define an array to hold edge information
edge_info = EdgeInfo[]

# Modified loop to store edge information
for i = 2:size(edge_data)[1]
    curr_edge = edge_data[i]
    curr_edge_dat = split(curr_edge, ',')
    edge_pieces = parse.(Int, curr_edge_dat[1:4])

    node_one = find_node_ind(edge_pieces[1], edge_pieces[2], node_data)
    node_two = find_node_ind(edge_pieces[3], edge_pieces[4], node_data)

    edge_length = parse(Float64, curr_edge_dat[5])

    broken_up = split(radii_data[i-1],",")
    radii_data_curr = parse.(Float64,broken_up)
    avg_radii = mean(radii_data_curr) * 1.83 * 10^(-6)

    push!(edge_info, EdgeInfo(i - 1, node_one, node_two, edge_length,
    (edge_pieces[1], edge_pieces[2]), 
    (edge_pieces[3], edge_pieces[4]), avg_radii))

end


cont = true
while cont
    node_degrees = Dict{Int, Int}()
    for edge in edge_info
        node_degrees[edge.in_node_idx] = get(node_degrees, edge.in_node_idx, 0) + 1
        node_degrees[edge.out_node_idx] = get(node_degrees, edge.out_node_idx, 0) + 1
    end
    nodes_to_remove = [node for (node, degree) in node_degrees if degree == 1]
    if isempty(nodes_to_remove) cont = false end
    edge_info = filter(edge -> !(edge.in_node_idx in nodes_to_remove || edge.out_node_idx in nodes_to_remove), edge_info)
end

#fix indexing
node_indices_in = [edge.in_node_idx for edge in edge_info]
node_indices_out = [edge.out_node_idx for edge in edge_info]

# Concatenate the two lists and then remove duplicates
node_indices = unique(vcat(node_indices_in, node_indices_out))
sort!(node_indices)

edge_indices = unique([edge.edge_idx for edge in edge_info])
sort!(edge_indices)

# Step 3: Map old indices to new indices
node_index_map = Dict(old_idx => new_idx for (new_idx, old_idx) in enumerate(node_indices))
edge_index_map = Dict(old_idx => new_idx for (new_idx, old_idx) in enumerate(edge_indices))

updated_edge_info = []

# Iterate through the existing edge_info_with_pos array
for edge in edge_info
    new_in_node_idx = node_index_map[edge.in_node_idx]
    new_out_node_idx = node_index_map[edge.out_node_idx]
    new_edge_idx = edge_index_map[edge.edge_idx]

    # Create a new EdgeInfo instance with the updated indices
    updated_edge = EdgeInfo(new_edge_idx, new_in_node_idx, new_out_node_idx, edge.edge_length, edge.in_node_pos, edge.out_node_pos, edge.avg_radii)

    # Add the updated edge to the new array
    push!(updated_edge_info, updated_edge)
end

edge_info = updated_edge_info

g = SimpleGraph(size(node_data)[1])

# Arrays to hold pixel length orders
pix_len_order = zeros(length(edge_info))
pix_len_order2 = zeros(length(edge_info))

# Array to hold edge nodes
edge_nodes = []
for i = 1:length(edge_info) push!(edge_nodes,[0;0]) end 
avg_radii = zeros(length(edge_info))

# Loop through the edge_info structure to add edges and calculate pixel lengths
for (i, edge) in enumerate(edge_info)
    node_one = edge.in_node_idx
    node_two = edge.out_node_idx
    edge_nodes[i] = [node_one; node_two]
    avg_radii[i] = edge.avg_radii

    add_edge!(g, node_one, node_two)

    # Assuming edge_length is already calculated and stored in edge_info
    pix_len_order[i] = edge.edge_length

    # Calculating the Euclidean distance between the nodes
    x1, y1 = edge.in_node_pos
    x2, y2 = edge.out_node_pos
    pix_len_order2[i] = ceil(sqrt((x1 - x2)^2 + (y1 - y2)^2))
end
