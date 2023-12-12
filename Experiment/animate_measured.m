%% Animate amplitudes on a Physarum network for measured values

filename = 'Sample one';
source = ['Data/' filename '/' filename];
destination = ['Movies/' filename '/' filename];

% Setup movie
writerObj = VideoWriter([destination '_movie_model.mp4'], 'MPEG-4');
open(writerObj);

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
colors = getPyPlot_cMap('RdYlBu', 500);



radii = readmatrix([source '_radii_time.txt']);

for i = 1:numEdges
    avg_radii(i) = mean(radii(i,:));
end

for i = 1:numEdges
    radii(i,:) = radii(i,:) / avg_radii(i);
end
for i = 1:numEdges
    for j = 1:size(radii,2)
        if radii(i,j) > 1.1
            radii(i,j) = 1.1;
        end
        if radii(i,j) < 0.9
            radii(i,j) = 0.9;
        end
    end
end

max_R = max(max(radii));
min_R = min(min(radii));


for X = 200:500%size(radii,2)
   disp(X)
    %Setup vessel colors
    map = zeros(numEdges,3);
    for j = 1:numEdges
        amp = radii(j,X);
        color_ind = round((amp-min_R)/(max_R-min_R) * 499)+1;
        map(j,:) = colors(color_ind,:);
    end

    curr_thick = mask;
%    for j = 1:Ne
%        curr = msk_edge_data(j).pix;
%        perc_change = radii(j,X)/avg_radii(j);
%        start_radius = avg_radii(j);
%        pixels_to_add = 2*ceil(start_radius * (perc_change-1));
%         if pixels_to_add > 0
%             if ~isempty(curr)
%                 for x = 1:length(curr)
%                     currx = curr(x,1);
%                     curry = curr(x,2);
%                     curr_thick(currx-pixels_to_add:currx+pixels_to_add,curry-pixels_to_add:curry+pixels_to_add) = 1;
%                 end
%             end
%         end
%         if pixels_to_add < 0
%             if ~isempty(curr)
%                 pixels_to_add = abs(pixels_to_add);
%                 for x = 1:length(curr)
%                     currx = curr(x,1);
%                     curry = curr(x,2);
%                     curr_thick(currx-pixels_to_add:currx+pixels_to_add,curry-pixels_to_add:curry+pixels_to_add) = 0;
%                 end
%             end
%         end
%     end

   

    fig_lab = figure(1); clf
    h = imagesc(diagram);
    alpha_data = logical(curr_thick);
    set(h, 'AlphaData', alpha_data)
    set(gca,'position',[0 0 1 1],'units','normalized')
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    axis equal
    colormap(map);
    
    drawnow
    writeVideo(writerObj, getframe(fig_lab));


end
close(writerObj);

