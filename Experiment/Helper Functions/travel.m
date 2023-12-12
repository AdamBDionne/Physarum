%function used to travel from a node 
function [graph, traveled] = travel(curr_node_x, curr_node_y, x, y, traveled, skel, graph, nodeList)
    x_end = 0; y_end = 0;
    cont = true;
    l = 1;
    edge_data(l,1) = x; edge_data(l,2) = y;

    while cont
        cont = false;
        keep = true;
        if l == 1
            for i = -1:1
                for j = -1:1
                    if i ~= 0 || j ~= 0
                        if skel(curr_node_x+i,curr_node_y+j) == 1 && traveled(curr_node_x+i,curr_node_y+j) == 0
                            if abs(curr_node_x+i-x) > 1 || abs(curr_node_y+j-y) > 1
                                traveled(curr_node_x,curr_node_y) = 1;
                                curr_node_x = curr_node_x+i; curr_node_y = curr_node_y+j;
                                cont = true;
                                l = l + 1;
                                edge_data(l,1) = curr_node_x; edge_data(l,2) = curr_node_y;
                                keep = false;
                                break
                            end
                        end
                    end
                end
                if ~keep
                    break
                end
            end
        end
        %check if next to a node and if so deal with this edge case
        if l > 2
            for i = -1:1
                for j = -1:1
                    if i~= 0 || j ~= 0 
                        if ismember([curr_node_x+i curr_node_y+j],nodeList,'rows')
                             traveled(curr_node_x,curr_node_y) = 1;
                             x_end = curr_node_x+i;
                             y_end = curr_node_y+j;
                             l = l + 1;
                            edge_data(l,1) = curr_node_x; edge_data(l,2) = curr_node_y;
                            keep = false;
                            break
                        end
                    end
                end
                if ~keep
                    break
                end
            end
        end
        if keep
            for i = -1:1
                for j = -1:1
                    if i ~= 0 || j ~= 0
                        if skel(curr_node_x+i,curr_node_y+j) == 1 && traveled(curr_node_x+i,curr_node_y+j) == 0
                            traveled(curr_node_x,curr_node_y) = 1;
                            curr_node_x = curr_node_x+i; curr_node_y = curr_node_y+j;
                            cont = true;
                            l = l + 1;
                            edge_data(l,1) = curr_node_x; edge_data(l,2) = curr_node_y;
                            break
                        end
                    end
                end
                if cont
                    break
                end
            end
        end
    end
    
    if x_end+y_end == 0
        for i = -2:2
            for j = -2:2
                if i~= 0 || j ~= 0 
                    if ismember([curr_node_x+i curr_node_y+j],nodeList,'rows')
                         traveled(curr_node_x,curr_node_y) = 1;
                         x_end = curr_node_x+i;
                         y_end = curr_node_y+j;
                         l = l + 1;
                            edge_data(l,1) = curr_node_x; edge_data(l,2) = curr_node_y;
                        break
                    end
                end
            end
        end
    end

    if x_end+y_end > 0
        curr_ind = length(graph) + 1;
        graph(curr_ind).node_one = [x y];
        graph(curr_ind).node_two = [x_end y_end];
        graph(curr_ind).length = l;
        graph(curr_ind).edge_data = edge_data;
    end
end
