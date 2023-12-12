%set params for loading in...
set(0,'DefaultFigureWindowStyle','docked')

filename = 'sample_two';
source = ['Raw data/' filename '/' filename];
destination = ['Processed data/' filename '/' filename];
data = [source '.tif'];


%load mask and graph
mask = imread([source '_mask_skel.png']);
filtered = medfilt2(mask(:,:,1), [9 9]);
filtered_again = imgaussfilt(filtered, 5);
BW = imbinarize(filtered_again,0.1);
mask = uint16(~BW);
mask = double(mask);


dims = size(mask);
graph = load([source '_graph.mat']);
graph = graph.graph;

%make voronoi cells
numEdges = length(graph);
%generate edge diagrams
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

clear edge_len
for j = 1:Ne
    node_one = graph(j).node_one;
    node_two = graph(j).node_two;
    edge_len(j) = sqrt((node_one(1) - node_two(1))^2 + (node_one(2) - node_two(2))^2);
end


y_min = 900; y_max = 3500; x_min = 1; x_max = 4600;

point_list_shifted(:,1) = point_list(:,1) - y_min + 1;
point_list_shifted(:,2) = point_list(:,2);

mask(:,4485:4600) = 0;
diagram(:,4486:4600) = 0;

if ~exist('msk_edge_data','var')
    for j = 1:Ne
        msk_edge_data(j).pix = [];
    end
    for i = 30:dims(1)-30
        for j = 30:dims(2)-30
            if (mask(i,j) ~= 0) && (sum(sum(mask(i-1:i+1,j-1:j+1))) ~= 9)
                crr = msk_edge_data(diagram(i,j)).pix;
                crr(end+1,1) = i;
                crr(end,2) = j;
                msk_edge_data(diagram(i,j)).pix = crr;
%                 edge_data(counter,1) = i;
%                 edge_data(counter,2) = j;
%                 edge_data(counter,3) = diagram(i,j);
%                 counter = counter+1;
            end
        end
    end
end

%% Color Map
    %map = [jet(128); 0 0 0];
    %caxis([0 pi+0.05])
    
%load processed phase data
Ne = 574;

%solution
sol_pts_dummy= readmatrix([source '_ODE_SOL.txt']);

%flows
flows_dummy = readmatrix([source '_FLOWS.txt']);

%translate 
translator = readmatrix([source '_translator.txt']);

flows = 0 * flows_dummy;
sol_pts = 0 * sol_pts_dummy;
for i = 1:Ne
    new_ind = floor(translator(i));
    sol_pts(:,new_ind) = sol_pts_dummy(:,i);
    sol_pts(:,Ne+new_ind) = sol_pts_dummy(:,Ne+i);
    flows(:,new_ind) = flows_dummy(:,i);
end
flows = transpose(flows); sol_pts = transpose(sol_pts);

colors = summer(256);
max_conc = max(max(sol_pts(Ne+1:end,:)));
min_conc = min(min(sol_pts(Ne+1:end,:)));
conc_norm = max([max_conc, abs(min_conc)]);

max_vol = max(max(sol_pts(1:Ne,:)));
min_vol = min(min(sol_pts(1:Ne,:)));
vol_norm = max([max_vol, abs(min_vol)]);

max_flow = max(max(flows));
min_flow = min(min(flows));
flow_norm = max([max_flow, abs(min_flow)]);


writerObj = VideoWriter([destination '_model_soln.mp4'], 'MPEG-4');
open(writerObj);

headWidth = 10;
headLength = 10;
LineLength = 1.0;



for t = 1:900
    disp(t)

    curr_vol = sol_pts(1:Ne,t);
    curr_conc = sol_pts(Ne+1:end,t);
    curr_flows = flows(:,t);

    %color edges
    map = zeros(numEdges+1,3);
    map(1,:) = [0 0 0];
    for j = 1:numEdges
        ind = ceil((curr_conc(j) - min_conc)/(max_conc-min_conc)*256);
        ind = max([ind,1]);
        map(j+1,:) = colors(ind,:);
        %curr_frame = curr_frame + ((~logical(curr_diagram - j)) .* modes(j,Ne-1));
    end

    %add thickness 
    curr_thick = mask;

    for j = 1:Ne
        curr = msk_edge_data(j).pix;
        num_to_add = floor((curr_vol(j)-min_vol)/(max_vol-min_vol)*15);
        if ~isempty(curr)
            for x = 1:length(curr)
                currx = curr(x,1);
                curry = curr(x,2);
                curr_thick(currx-num_to_add:currx+num_to_add,curry-num_to_add:curry+num_to_add) = 1;
            end
        end
    end
    
    %quiver
    %find arrows starting position
    clear arrow_x_pos; clear arrow_y_pos;
    clear arrow_direc_x; clear arrow_direx_y;
    for j = 1:Ne
        node_one = graph(j).node_one;
        node_two = graph(j).node_two;
        dx = node_one(1) - node_two(1);
        dy = node_one(2) - node_two(2);
        theta = atan(dy/dx);
       norm = sqrt((- node_one(1) + node_two(1))^2 + (- node_one(2) + node_two(2))^2);
        
        edge_flow = curr_flows(j);
        normed_flow = edge_flow / flow_norm;
        if abs(normed_flow) < 0.6
            normed_flow = sign(normed_flow)*0.6;
        end
  
        arrow_len = abs(normed_flow * sqrt((node_one(1)-node_two(1))^2 + (node_two(2) - node_one(2))^2));
        if edge_flow > 0
            %flow from node two to node one
            %shift towards node two
            arrow_direc_x(j) = (node_one(1) - node_two(1))*abs(arrow_len)/norm;
            arrow_direc_y(j) = (node_one(2) - node_two(2))*abs(arrow_len)/norm;
                if node_one(1) > node_two(1)
                    arrow_x_pos(j) = point_list_shifted(j,1) - abs(cos(theta))*arrow_len/2;
                else
                    arrow_x_pos(j) = point_list_shifted(j,1) + abs(cos(theta))*arrow_len/2;
                end
                
                if node_one(2) > node_two(2)
                    arrow_y_pos(j) = point_list_shifted(j,2) - abs(sin(theta))*arrow_len/2;
                else
                    arrow_y_pos(j) = point_list_shifted(j,2) + abs(sin(theta))*arrow_len/2;
                end
        else
            %flow from node one to node two
            %shift towards node one
         
            arrow_direc_x(j) = (- node_one(1) + node_two(1))*abs(arrow_len)/norm;
            arrow_direc_y(j) = (- node_one(2) + node_two(2))*abs(arrow_len)/norm;
                if node_one(1) > node_two(1)
                    arrow_x_pos(j) = point_list_shifted(j,1) + abs(cos(theta))*arrow_len/2;
                else
                    arrow_x_pos(j) = point_list_shifted(j,1) - abs(cos(theta))*arrow_len/2;
                end
                
                if node_one(2) > node_two(2)
                    arrow_y_pos(j) = point_list_shifted(j,2) + abs(sin(theta))*arrow_len/2;
                else
                    arrow_y_pos(j) = point_list_shifted(j,2) - abs(sin(theta))*arrow_len/2;
                end
        end
        
    end
    
    figure(2); imagesc(mask); 
    set(gca,'XDir','normal')
    set(gca,'YDir','normal')
    axis equal
    hq = quiver(arrow_y_pos,arrow_x_pos,arrow_direc_y,arrow_direc_x);
    figure(1); clf; ax = axes();
    ax.Colormap = map;
    h = imagesc(diagram(y_min:y_max,x_min:x_max));
    alpha_data = logical(curr_thick(y_min:y_max,x_min:x_max));
    set(h, 'AlphaData', alpha_data)
    set(gca,'position',[0 0 1 1],'units','normalized')
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'XDir','normal')
    set(gca,'YDir','normal')
    axis equal
    hold on
    ax2 = axes('Visible','off','HandleVisibility','off'); 
    cmap = colormap(map);
    ax2.Colormap = cmap;
    set(gca,'YDir','normal')
    set(gca,'XDir','normal')
    axis equal
    
   
    %get the data from regular quiver
    U = hq.UData;
    V = hq.VData;
    X = hq.XData;
    Y = hq.YData;

    for ij = 1:length(X)
        headWidth = 5;
        ah = annotation('arrow',...
            'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth);
        set(ah,'parent',gca);
        set(ah,'position',[X(1,ij) Y(1,ij) LineLength*U(1,ij) LineLength*V(1,ij)]);
    end
    
    drawnow
    writeVideo(writerObj, getframe(gcf));
end

close(writerObj);

