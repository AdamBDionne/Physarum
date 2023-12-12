

function [inside] = circle_inside(pix_x, pix_y, I, radius, thresh, circ_pts_library)
    inside = true;
    dims = size(I);
    circ_pts = circ_pts_library(radius).pts;
    circ_pts(:,1) = circ_pts(:,1) + pix_x;
    circ_pts(:,2) = circ_pts(:,2) + pix_y;
    for j = 1:length(circ_pts)
        x = circ_pts(j,1);
        y = circ_pts(j,2);
        if x > 0 && y > 0 && x <= dims(1) && y <= dims(2) && I(x,y) > thresh
            inside = false;
            break
        end
    end
end