%% Overlays particles form nutrient dispersal simulation

%Load in data
amp_t = load("amp_t.txt");
part_pos_t = load("fin_part_pos.mat");
part_pos_t = part_pos_t.fin_part_pos;
%part_pos_t(:,1,:) = 2*part_pos_t(:,1,:);

colors = getPyPlot_cMap('RdYlBu', 500, [], '"C:\\Users\\abd12\\AppData\\Local\\Microsoft\\WindowsApps\\PythonSoftwareFoundation.Python.3.10_qbz5n2kfra8p0\\python.exe"');

%%% load in data for network
filename = 'sample two';
source = ['Data/' filename '/' filename];

graph = load([source '_graph.mat']);
graph = graph.graph;

crop = load([source '_crop.txt']);
y_min = crop(1); y_max = crop(2);
x_min = crop(3); x_max = crop(4);
numEdges = length(graph); Ne = numEdges;

translator = readmatrix([source, '_translator.txt']);

% Mask
mask = imread([source '_mask.png']);
filtered = medfilt2(mask(:,:,1), [20 20]);
filtered_again = imgaussfilt(filtered, 5);
BW = imbinarize(filtered_again,0.1);
mask = uint16(~BW);
mask = double(mask);

mask = mask(y_min:y_max,x_min:x_max);
dims = size(mask);

% Make voronoi cells for coloring
if ~exist('diagram','var')
    point_list = zeros(numEdges,2);

    for i = 1:numEdges
        edgeData = graph(i).edge_data;
        point_list(i,1) = mean(edgeData(:,1));
        point_list(i,2) = mean(edgeData(:,2));
    end

    diagram = zeros(dims(1), dims(2));

    for i = 1:dims(1)
        for j = 1:dims(2)
            [~,ind] = min((point_list(:,1) - i).^2 + (point_list(:,2) - j).^2);
            diagram(i,j) = ind;
        end
    end
end


writerObj = VideoWriter('Nutrient Transport/results/movie_transport.mp4', 'MPEG-4');
writerObj.FrameRate = 10;
open(writerObj);
for frame = 1:200
    phase = amp_t(frame,:);
    part_pos = part_pos_t(:,:,frame);

    %translate 
    temp = phase;
    temp2 = part_pos;
    for i = 1:numEdges
        phase(abs(translator(i))) = temp(i);
        if size(part_pos,2) > 1
            part_pos(abs(translator(i)),1) = temp2(i,1);
            part_pos(abs(translator(i)),2) = abs(translator(temp2(i,2)));
        else
            part_pos(abs(translator(i))) = temp2(i);
        end
    end
    
    % Color Map
    max_amp = max(phase);
    min_amp = min(phase);
    
    %Setup phase colors
    map = zeros(numEdges,3);
    %map(1,:) = [0 0 0];
    for j = 1:numEdges
        amp = phase(j);
        color_ind = round((amp-min_amp)/(max_amp-min_amp) * 499)+1;
        map(j,:) = colors(color_ind,:);
    end
    
    fig_lab = figure(1); clf
    h = imagesc(diagram);
    alpha_data = logical(mask);
    set(h, 'AlphaData', alpha_data)
    set(gca,'position',[0 0 1 1],'units','normalized')
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    axis equal
    colormap(map);
    hold all
    
    
    %plot particles
    
    %max + min y pos
    maxy = 0;
    miny = 1e6;
    for i = 1:Ne
        node_one = graph(i).node_one;
        node_two = graph(i).node_two;
        mid_y = (node_one(1) + node_two(1))/2;
        if mid_y > maxy
            maxy = mid_y;
        end
        if mid_y < miny
            miny = mid_y;
        end
    end
    
    radius = 7;
    outline_width = 2;
    turbo_map = turbo(64);
    for i = 1:Ne
        node_one = graph(i).node_one;
        node_two = graph(i).node_two;
        pos = part_pos(i,1);
        y_mid = (node_one(1) + node_two(1))/2;
        % Calculate the particle's position based on part_pos
        if size(part_pos,2) > 1
            node_one = graph(part_pos(i,2)).node_one;
            node_two = graph(part_pos(i,2)).node_two;
            particle_position = node_one + pos * (node_two - node_one);
        else
            particle_position = node_one + pos * (node_two - node_one);
        end
    
        % Calculate the color based on the Y position of (node_one + node_two)/2
       
        color = turbo_map(floor(1+63*(y_mid-miny)/(maxy-miny)),:);
        scatter(particle_position(2), particle_position(1), (radius + outline_width)^2, 'k', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'none');
        scatter(particle_position(2), particle_position(1), radius^2, color, 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'none');
    end
    drawnow
    writeVideo(writerObj, getframe(fig_lab));
end

close(writerObj);