%given skeleton matrix, find the largest edge
function [edge_pixs] = find_largest_edge(skel)
    [graph, ~] = construct_network(skel); 
    max = 0;
    for i = 1:length(graph)
        if graph(i).length > max
            max = graph(i).length;
            ind = i;
        end
    end
    edge_pixs = graph(ind).edge_data;
end

