%takes a skeleton and returns a cell array with graph data
function [graph, nodeList] = construct_network(skel)
    %find nodes from skeleton based on adjacency
    c =1;
    dims = size(skel);
    traveled = 0.*skel;
    graph = struct([]);

    for i = 2:dims(1)-1
        for j= 2:dims(2)-1
            if skel(i,j) == 1 && (sum(sum(skel(i-1:i+1,j-1:j+1))) > 3 || sum(sum(skel(i-1:i+1,j-1:j+1))) == 2) && sum(sum(traveled(i-1:i+1,j-1:j+1))) == 0
              nodeList(c,1) = i; nodeList(c,2) = j;
              traveled(i,j) = 1;
              c = c+1;
            end
        end
    end


    %go through each node, and then travel from that node
    %until we reach another node. this makes an edge, and we mark the 
    %edge and nodes as travelled
    for k = 1:length(nodeList)
        curr_node_x = nodeList(k,1);
        curr_node_y = nodeList(k,2);
        for i = -1:1
            for j = -1:1
                if i ~= 0 || j ~= 0
                    if skel(curr_node_x+i,curr_node_y+j) == 1 && traveled(curr_node_x+i,curr_node_y+j) == 0
                        [graph, traveled] = travel(curr_node_x+i, curr_node_y+j, curr_node_x, curr_node_y, traveled, skel, graph, nodeList);
                    end
                end
            end
        end
    end
end

