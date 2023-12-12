%% Measure the radii of each vessel in network
%% and animate microscopy movie

filename = 'sample_spring_seven';
source = ['Raw data/' filename '/' filename];
destination = ['Processed data/' filename '/' filename];
data = [source '.tif'];

info = imfinfo(data);
numFrames = length(info);
start = 70;

fps = 5;

%%% load in data

graph = load([source '_graph.mat']);
graph = graph.graph;

crop = load([source '_crop.txt']);
y_min = crop(1); y_max = crop(2);
x_min = crop(3); x_max = crop(4);

%%% measure radii & animate
numEdges = length(graph);
map = colormap(gray);

% Setup movie
writerObj = VideoWriter([destination '_movie_s.mp4'], 'MPEG-4');
open(writerObj);

radii_results = zeros(numFrames-start+1, numEdges);

%Generate frames
for ind = start:numFrames
    tic
    disp(['----- Frame ' int2str(ind) ' -----'])

    raw = imread(data,'Index',ind);
    raw = raw(y_min:y_max,x_min:x_max,1);

    t = imgaussfilt(raw,[50 50]);
    t2 = (raw-t) == 0; 
    t3 = medfilt2(t2,[20 20]);
    t3 = ~t3;

    for i = 1:numEdges
        edge_pixs = graph(i).edge_data;
        if ind == start
            radius = measure_radius(edge_pixs, t3, 1, 0.9, circ_pts_library);
        else
            radius = measure_radius(edge_pixs, t3, floor(radii_results(ind-start+1,i)), 0.9, circ_pts_library);
        end
        radii_results(ind-start+1,i) = radius;
        if isnan(radius)
            disp(i)
        end

    end

    raw(end:end+300,:) = max(max(raw)); id = ind-start+1; 
    fig_lab = figure(1); clf; imagesc(raw); colormap(map); axis equal;

    set(gca,'position',[0 0 1 1],'units','normalized')
    set(gca,'XTick',[])
    set(gca,'YTick',[])

    % Time bar
    mins = floor((id*fps)/60);
    if mins == 1
        curr_time = [int2str(mins) ' minute'];
    else
        curr_time = [int2str(mins) ' minutes'];
    end
    annotation(fig_lab,'textbox',...
    [0.7 0.07 0.3 0.05],...
    'String',curr_time,...
    'FontSize',40,...
    'FitBoxToText','off',...
    'EdgeColor','none');

    % Seconds
    secs = mod(id,60/fps)*fps;

     annotation(fig_lab,'textbox',...
    [0.7 0.035 0.3 0.05],...
    'String',[int2str(secs) ' seconds'],...
    'FontSize',40,...
    'FitBoxToText','off',...
    'EdgeColor','none');

    % Scale bar
    xlen = x_max-x_min;
    scale_bar_len = (1000/1.83) / (xlen);

    annotation(fig_lab,'rectangle',...
    [0.05 0.05 scale_bar_len 0.025],...
    'FaceColor',[0 0 0]);

    x_label_pos = (1/2)*(0.1+scale_bar_len) - 0.1;

    annotation(fig_lab,'textbox',...
    [x_label_pos 0.07 0.25 0.05],...
    'String','1 millimeter',...
    'FontSize',40,...
    'EdgeColor','none');

    % Save
    drawnow
    writeVideo(writerObj, getframe(fig_lab));
    toc
end

% Save results
close(writerObj);
save([source '_radii.mat'],'radii_results')

