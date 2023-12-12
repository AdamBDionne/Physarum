%Still frame of phases


savefilename = "real_largeproj.png";
loadinname = "frame_modelC.txt";
colors = getPyPlot_cMap('RdYlBu', 500, [], '"C:\\Users\\abd12\\AppData\\Local\\Microsoft\\WindowsApps\\PythonSoftwareFoundation.Python.3.10_qbz5n2kfra8p0\\python.exe"');
%colors = viridis(500);

filename = 'sample_spring_seven';
source = ['Raw data/' filename '/' filename];
destination = ['Processed data/' filename '/' filename];

% local_order = load([source '_local_order.txt']);
% global_order = load([source '_global_order.txt']);
% [~,time] = max(local_order(100:end-100) + (1 - global_order(100:end-100)));
% %[~,time] = min(global_order(100:end-100));
% time = time + 100;
% time = 302;

%%% load in data
graph = load([source '_graph.mat']);
graph = graph.graph;

crop = load([source '_crop.txt']);
y_min = crop(1); y_max = crop(2);
x_min = crop(3); x_max = crop(4);

% %compute phases
% radii_signal = load([source '_radii_time.txt']);
% 
 numEdges = length(graph);
% phase_signal = 0*radii_signal;
% 
% for j = 1:numEdges
%     phase_signal(j,:) = angle(hilbert(radii_signal(j,:) - mean(radii_signal(j,:))));
% end

%phase = phase_signal(:,time);

%phase = load(['Processed Data/' filename '/' filename '_MODES.txt']);
%phase = phase(:,numEdges-0);

phase = load(loadinname);
% for j =1:size(phase,1)
%     if phase(j) < -0.4
%         phase(j) = -0.4;
%     end
%     if phase(j) > 0.4
%         phase(j) = 0.4;
%     end
% end

translator = readmatrix([source, '_translator.txt']);
%mode = mode(abs(translator));
temp = phase;
for i = 1:numEdges
    phase(abs(translator(i))) = temp(i);
end

% t = sort(phase);
% newmax = t(end-10); newmin = t(10);
% for ind = 1:length(phase)
%     if phase(ind) > newmax
%         phase(ind) = newmax;
%     end
%     if phase(ind) < newmin
%         phase(ind) = newmin;
%     end
% end

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


% Color Map


%diagram_cropped = diagram(y_min:y_max,x_min:x_max);
diagram_cropped = diagram;
max_amp = max(phase);
min_amp = min(phase);

%Setup phase colors
map = zeros(numEdges,3);
%map(1,:) = [0 0 0];
for j = 1:numEdges
    amp = phase(j);
    color_ind = round( (amp-min_amp)/(max_amp-min_amp) * 499)+1;
    map(j,:) = colors(color_ind,:);
end

fig_lab = figure(); clf
h = imagesc(diagram_cropped);
alpha_data = logical(mask);
set(h, 'AlphaData', alpha_data)
set(gca,'position',[0 0 1 1],'units','normalized')
set(gca,'XTick',[])
set(gca,'YTick',[])
axis equal
colormap(map);


%saveas(gcf,[destination 'phase_global.png'])
%dlmwrite([source '_phase_global_time.txt'],time)

saveas(gcf,savefilename)



