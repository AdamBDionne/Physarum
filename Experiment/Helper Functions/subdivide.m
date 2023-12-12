function [graph, nodeList] = subdivide(graph, nodeList)
    
    %generate circle points
    if ~exist('circ_pts_library','var')
        circ_pts_library = struct([]);
        for r = 1:200
            counter = 1;
            circ_pts = [r 0];
            for theta = [0 (pi/2) 0:(pi/(10*r)):pi/2]
                x = sign(cos(theta))*floor(abs(r*cos(theta)));
                y = sign(sin(theta))*floor(abs(r*sin(theta)));
                new_pt = [x y];
                if  ~sum(ismember(circ_pts,new_pt,'rows'))
                    circ_pts(counter,:) = new_pt;
                    counter = counter+1;
                    if new_pt(1) ~= 0
                        circ_pts(counter,:) = [-new_pt(1) new_pt(2)];
                        counter = counter + 1;
                    end
                    if new_pt(2) ~= 0
                        circ_pts(counter,:) = [new_pt(1) -new_pt(2)];
                        counter = counter + 1;
                    end
                    if (new_pt(1) ~= 0) && (new_pt(2) ~= 0)
                        circ_pts(counter,:) = [-new_pt(1) -new_pt(2)];
                        counter = counter + 1;
                    end
                end
            end

            circ_pts_library(r).pts = circ_pts;
        end
    end
    
     
    %subdivide to make sure edges are accurate to shape

    stack = java.util.Stack();
    for i = 1:length(graph)
        stack.push(i);
    end

    while ~stack.empty()
        i = stack.pop();
        node_one = graph(i).node_one;
        node_two = graph(i).node_two;
        edge_data = graph(i).edge_data;
        len = graph(i).length;
        
        
        %make line connecting two nodes 
        if node_one(1) ~= node_two(1)
            acoeff = ( node_one(2) - node_two(2) )/ (node_two(1) - node_one(1));
            bcoeff = 1;
            ccoeff = -(node_one(2)+acoeff*node_one(1));
        else
            acoeff = 1;
            bcoeff = (node_one(1) - node_two(1))/(node_two(2) - node_one(2));
            ccoeff = -(node_one(1)+bcoeff*node_one(2));
        end
        
        norm = sqrt(acoeff^2+bcoeff^2);
        [val,ind] = max(abs(acoeff*edge_data(:,1) + bcoeff*edge_data(:,2)+ccoeff)/norm);
        
        
        
        if val > 20 || len > 150
            %now we will split this edge
            if ind == 1 || ind == len
                ind = floor(len / 2);
            end

            node_three = edge_data(ind,:);
            if sum(node_one == edge_data(1,:)) == 2
                edge_data_one = edge_data(1:ind,:);
                edge_data_two = edge_data(ind:end,:);
            else
                edge_data_two = edge_data(1:ind,:);
                edge_data_one = edge_data(ind:end,:);
            end

            graph(i).node_one = node_one;
            graph(i).node_two = node_three;
            graph(i).edge_data = edge_data_one;
            graph(i).length = length(edge_data_one);

            curr_ind = length(graph) + 1;
            graph(curr_ind).node_one = node_two;
            graph(curr_ind).node_two = node_three;
            graph(curr_ind).edge_data = edge_data_two;
            graph(curr_ind).length = length(edge_data_two);
            stack.push(i);
            stack.push(curr_ind);

            nodeList(length(nodeList)+1,:) = node_three;
        end
    end


    %find volume of each edge 
%         raw= imread(data,'Index',1);
%         raw = raw(y_min:y_max,x_min:x_max,1);
%         raw = imresize(raw,0.5);
% 
%         %process 
%         I_smooth = rog_smooth(raw, 0.003, 3, 3, 3);
%         I = imbinarize(I_smooth);
%          edge_vols = zeros(length(graph),1);
%     
%     for i = 1:length(graph)
%             edge_pixs = graph(i).edge_data;
%             radius = measure_radius(edge_pixs, BW, 1, 0.9, circ_pts_library);
%             leng = graph(i).length;
%             vol = pi*radius^2*leng;
%             graph(i).vol = vol;
%             edge_vols(i) = vol;
%             hold all; plot(mean(edge_pixs(:,2)), mean(edge_pixs(:,1)),'*r');
%             viscircles([mean(edge_pixs(:,2)), mean(edge_pixs(:,1))], radius);
%     end
%     
%     med_vol = median(edge_vols);
%     
%     stnd_vol = 2*med_vol;
%     
%     new_graph = struct();
%     counter = 1;
%     for i = 1:length(graph)
%         curr_vol = graph(i).vol;
%         num_segs = round(curr_vol/stnd_vol);
%         % extend / reduce to stnd vol
%         % one edge
%         if num_segs < 2
%             new_graph(counter).node_one = graph(i).node_one;
%             new_graph(counter).node_two = graph(i).node_two;
%             new_graph(counter).edge_data = graph(i).edge_data;
%             new_graph(counter).length = stnd_vol / (graph(i).vol / graph(i).length);
%             new_graph(counter).vol = stnd_vol;
%             counter = counter+1;
%         % multiple edges
%         else
%             new_nodes = [1 floor((1:num_segs-1)*graph(i).length/num_segs) graph(i).length];
%             edge_data = graph(i).edge_data;
%             for j = 1:num_segs
%                 new_graph(counter).node_one = edge_data(new_nodes(j),:);
%                 new_graph(counter).node_two = edge_data(new_nodes(j+1),:);
%                 new_graph(counter).edge_data = edge_data(new_nodes(j):new_nodes(j+1),:);
%                 new_graph(counter).length = stnd_vol / (graph(i).vol / graph(i).length);
%                 new_graph(counter).vol = stnd_vol;
%                 if j ~= num_segs
%                      nodeList(length(nodeList)+1,:) = edge_data(new_nodes(j+1),:);
%                 end
%                 counter = counter + 1;
%             end
%         end
%     end
%     
%     graph = new_graph;
   
end