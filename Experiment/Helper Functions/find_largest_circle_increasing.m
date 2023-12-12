
function [radius] = find_largest_circle_increasing(pix_x, pix_y, I, radius, thresh, circ_pts_library)
    while radius < 200
        if ~circle_inside(pix_x, pix_y, I, radius, thresh, circ_pts_library)
            break
        else
            radius = radius + 1;
        end
    end
    
    if radius >= 200
        radius = NaN;
        disp('big')
    end
end