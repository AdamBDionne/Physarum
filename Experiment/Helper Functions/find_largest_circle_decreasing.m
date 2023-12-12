
function [radius] = find_largest_circle_decreasing(pix_x, pix_y, I, radius, thresh, circ_pts_library)
    while radius > 1
        if circle_inside(pix_x, pix_y, I, radius, thresh, circ_pts_library)
            break
        else
            radius = radius - 1;
        end
    end
    
    if radius <= 1
        radius = NaN;
        %disp('small')
    end

end