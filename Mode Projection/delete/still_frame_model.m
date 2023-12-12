%Still frame of model

still_ind = 20;

filename = 'sample_spring_three';
source = ['Raw data/' filename '/' filename];
destination = ['Processed data/' filename '/' filename];

%%% load in data
graph = load([source '_graph.mat']);
graph = graph.graph;

crop = load([source '_crop.txt']);
y_min = crop(1); y_max = crop(2);
x_min = crop(3); x_max = crop(4);

%%% measure radii & animate
numEdges = length(graph);
Ne = numEdges;


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

if ~exist('msk_edge_data','var')
    for j = 1:numEdges
        msk_edge_data(j).pix = [];
    end
    for i = 30:dims(1)-30
        for j = 30:dims(2)-30
            if (mask(i,j) ~= 0) && (sum(sum(mask(i-1:i+1,j-1:j+1))) ~= 9)
                crr = msk_edge_data(diagram(i,j)).pix;
                crr(end+1,1) = i;
                crr(end,2) = j;
                msk_edge_data(diagram(i,j)).pix = crr;
            end
        end
    end
end


% Color Map
colors = viridis(500);
%colors = getPyPlot_cMap('viridis', 500);


% Load in solution
sol = readmatrix([source, '_ODE_SOL.txt']);
model_sol = sol(still_ind,:);

max_c = max(model_sol(numEdges+1:end))-0.05;
min_c = min(model_sol(numEdges+1:end));
max_R = max(model_sol(1:numEdges));

radii = readmatrix([source '_radii_time.txt']);
for i = 1:numEdges
    avg_radii(i) = mean(radii(i,:));
end


%translate from matlab indexing to julia indexing 
translator = readmatrix([source, '_translator.txt']);
%mode = mode(abs(translator));
temp = model_sol;
for i = 1:numEdges
    model_sol(abs(translator(i))) = temp(i);
    model_sol(numEdges+abs(translator(i))) = temp(numEdges+i);
end


%Setup concentration colors
map = zeros(numEdges,3);
for j = 1:numEdges
    conc = model_sol(numEdges+j);
    color_ind = round( (conc-min_c)/(max_c-min_c) * 499)+1;
    map(j,:) = colors(color_ind,:);
end

curr_thick = mask;
for j = 1:Ne
    curr = msk_edge_data(j).pix;
    perc_change = model_sol(j);
    start_radius = avg_radii(j);
    pixels_to_add = 2*ceil(start_radius * (perc_change-1));
    if pixels_to_add > 0
        if ~isempty(curr)
            for x = 1:length(curr)
                currx = curr(x,1);
                curry = curr(x,2);
                curr_thick(currx-pixels_to_add:currx+pixels_to_add,curry-pixels_to_add:curry+pixels_to_add) = 1;
            end
        end
    end
    if pixels_to_add < 0
        if ~isempty(curr)
            pixels_to_add = abs(pixels_to_add);
            for x = 1:length(curr)
                currx = curr(x,1);
                curry = curr(x,2);
                curr_thick(currx-pixels_to_add:currx+pixels_to_add,curry-pixels_to_add:curry+pixels_to_add) = 0;
            end
        end
    end
end

%Setup expansions and contractions 


fig_lab = figure(start_ind); clf
h = imagesc(diagram);
alpha_data = logical(curr_thick);
set(h, 'AlphaData', alpha_data)
set(gca,'position',[0 0 1 1],'units','normalized')
set(gca,'XTick',[])
set(gca,'YTick',[])
axis equal
colormap(map);


saveas(gcf,[destination '_model_still_' int2str(still_ind) '.png'])