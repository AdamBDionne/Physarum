
function [radius] = measure_radius(edge_pixs, I, previous_radius, thresh, circ_pts_library)

    num_pixs = size(edge_pixs,1);
    counter = 1;
    radii = [];
    for i = 1:num_pixs
        if previous_radius == 0 || isnan(previous_radius)
            radius = find_largest_circle_increasing(edge_pixs(i,1), edge_pixs(i,2), I, 1, thresh, circ_pts_library);
        else
            if circle_inside(edge_pixs(i,1), edge_pixs(i,2), I, previous_radius, thresh, circ_pts_library)
                radius = find_largest_circle_increasing(edge_pixs(i,1), edge_pixs(i,2), I, previous_radius+1, thresh, circ_pts_library);
            else
                radius = find_largest_circle_decreasing(edge_pixs(i,1), edge_pixs(i,2), I, previous_radius-1, thresh, circ_pts_library);
            end
        end

        if ~isnan(radius)
            radii(counter) = radius;
            counter = counter+1;
        end
    end
    
    filtered_radii = rmoutliers(radii);
%     filtered_radii = transpose(filtered_radii);
%     num_pts = length(filtered_radii);
%     mid_pt = (num_pts)/2;
%     weight = zeros(num_pts,1);
%     sigma = (num_pts/1.5)^2;
%     for indx = 1:num_pts
%         weight(indx) = exp(-(indx-mid_pt)^2/sigma); 
%     end
%     weight = weight ./ sum(weight);
%     radius = sum(weight .* filtered_radii);
    radius = mean(filtered_radii);
end