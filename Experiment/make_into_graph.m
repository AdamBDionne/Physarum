%Here we skeletonize a reference physarum image, and then decompose the
%skeleton into a graph. 
save_res = false;

%% MAKING SKELETON

%load in mask
mask = imread([source '_mask.png']);

% Crop
figure();
I = imread(data,'Index',1);
imagesc(mask);
disp('Select top left corner...')
[x_min, y_min] = ginput(1);

disp('Select bottom right corner...')
[x_max, y_max] = ginput(1);

x_min = floor(x_min); y_min = floor(y_min); 
x_max = floor(x_max); y_max = floor(y_max); 

mask = mask(y_min:y_max, x_min:x_max,:); 


%filter
filtered = medfilt2(mask(:,:,1), [9 9]);
filtered_again = imgaussfilt(filtered, 5);

%binarize
BW = imbinarize(filtered_again);

%make skeleton
skel = bwskel(~BW);

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
    

%% MAKING NETWORK 
   
[graph, nodeList] = construct_network(skel);

%see network first pass
figure(); imagesc(BW + skel); hold all; plot(nodeList(:,2),nodeList(:,1),'.','MarkerSize',30)
for i = 1:length(graph)
   
    node_one = graph(i).node_one;
    node_two = graph(i).node_two;
    edge_data = graph(i).edge_data;
    plot([node_one(2) node_two(2)],[node_one(1) node_two(1)],'LineWidth',5)
    
    %label each edge
    txt = int2str(i);
    text((1/2)*(node_one(2)+node_two(2)),(1/2)*(node_one(1)+node_two(1)),txt,'FontSize',20)
 
end

%remove extraneous edges by eye
to_remove = input("What edges do you want to remove? Format [2 4 5])");

for i = 1:length(to_remove)
    edge_data = graph(to_remove(i)).edge_data;
    for j = 2:length(edge_data)-1
        skel(edge_data(j,1), edge_data(j,2)) = 0;
    end
end

%remove the nodes that are now floating in skel..

CC = bwconncomp(skel);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);

skel = 0.*skel;
skel(CC.PixelIdxList{idx}) = 1;

%do it all again...
[graph, nodeList] = construct_network(skel);

%remove junction pixels
figure(); imagesc(BW);
for i = 1:length(graph)
    edge_data = graph(i).edge_data; 
    start_pix = graph(i).node_one;
    end_pix = graph(i).node_two;
    %find circle for both starting and ending pixel
    
    start_circle_rad = find_largest_circle_increasing(start_pix(1), start_pix(2), BW, 1, 0.9, circ_pts_library);
    end_circle_rad = find_largest_circle_increasing(end_pix(1), end_pix(2), BW, 1, 0.9, circ_pts_library);
    
    outside_junction_edge_data = [];
    for j = 1:size(edge_data,1)
        pix_x = edge_data(j,1);
        pix_y = edge_data(j,2);
        if sqrt((pix_x -start_pix(1))^2+(pix_y-start_pix(2))^2) > start_circle_rad
            if sqrt((pix_x-end_pix(1))^2+(pix_y-end_pix(2))^2) > end_circle_rad
                outside_junction_edge_data(end+1,:) = [pix_x pix_y];
            end
        end
    end
    if isempty(outside_junction_edge_data)
        outside_junction_edge_data = edge_data;
    end
    graph(i).edge_data = graph(i).edge_data;
    hold all; plot(outside_junction_edge_data(:,2),outside_junction_edge_data(:,1));
    viscircles([start_pix(2) start_pix(1)], start_circle_rad); viscircles([end_pix(2) end_pix(1)], end_circle_rad);
end


%Now, go through each edge and split up until edges have variation
%underneath some set threshhold

[graph, nodeList] = subdivide(graph, nodeList);

%remove edge case where an edge is from a node to itself
to_remove = [];
for i = 1:length(graph)
    node_one = graph(i).node_one;
    node_two = graph(i).node_two;
    if sum(node_one == node_two) == 2
        to_remove(end+1) = i;
    end
end
for i = 1:length(to_remove)
    graph(to_remove(i)-i+1) = [];
end

%see updated network
fig = figure(); clf; imagesc(BW); hold all %imagesc(skel); hold all; plot(nodeList(:,2),nodeList(:,1),'.','MarkerSize',30)
for i = 1:length(graph)
    node_one = graph(i).node_one;
    node_two = graph(i).node_two;
    edge_data = graph(i).edge_data;
    plot([node_one(2) node_two(2)],[node_one(1) node_two(1)],'LineWidth',5,'Color',[30/255 138/255 222/255])
    plot(node_one(2),node_one(1),'.','MarkerSize',35,'Color',[209/255 70/255 70/255])
    plot(node_two(2),node_two(1),'.','MarkerSize',35,'Color',[209/255 70/255 70/255])
    %pause(0.01)
    
    %label each edge
    txt = int2str(i);
    text((1/2)*(node_one(2)+node_two(2)),(1/2)*(node_one(1)+node_two(1)),txt,'FontSize',20)
end



%% save skeleton, network, and plot
network = struct([]);
for i = 1:length(graph)
    network(i).node_one = graph(i).node_one;
    network(i).node_two = graph(i).node_two;
    network(i).length = graph(i).length;
end

nodeList_cmpl = []; clear i;
for j = 1:length(graph)
    n1 = graph(j).node_one;
    n2 = graph(j).node_two; 
    cmpl1 = n1(1) + (n1(2)*i);
    cmpl2 = n2(1) + (n2(2)*i);
    if sum(nodeList_cmpl == cmpl1) == 0
        nodeList_cmpl(end+1) = cmpl1;
    end
    if sum(nodeList_cmpl == cmpl2) == 0
        nodeList_cmpl(end+1) = cmpl2;
    end
end

nodeList = [];
for j = 1:length(nodeList_cmpl)
    nodeList(j,1) = real(nodeList_cmpl(j));
    nodeList(j,2) = imag(nodeList_cmpl(j));
end


if save_res
    saveas(gcf,[destination '_overlay_blank.png'])
    writetable(struct2table(network),[source '_network_edges.csv'])
    dlmwrite([source '_network_nodes.txt'], nodeList, 'delimiter','\t')
    dlmwrite([source '_skeleton.txt'], skel, 'delimiter','\t')
    save([source '_graph.mat'],'graph')
    dlmwrite([source '_crop.txt'],[y_min;y_max;x_min;x_max],'delimiter','\t')
end





