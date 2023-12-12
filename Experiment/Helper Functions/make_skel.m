%make skeleton from Physarum image
function [skel] = make_skel(I, mask, dims)


    bgillum = imgaussfilt(I,100);
    I_in = I - bgillum;
    BW = ~imbinarize(I_in,0) .* logical(mask);
    
%     imagesc(I_in .- (0.*I_in));
%     %binarize
%     BW = imbinarize(I,0.2);
%     BW = BW(:,:,1);
%     BW = ~BW;
%     BW = BW .* logical(mask);

    %filter noise
    filt = medfilt2(BW, [15 15]);

    %take largest connected comp.
    tubes = zeros(dims(1),dims(2));
    CC = bwconncomp(filt);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    tubes(CC.PixelIdxList{idx}) = 1;

    %remove small holes
    tubes = ~bwareaopen(~tubes, 500);

    %smooth
    smoothed = imgaussfilt(double(tubes),10);

    %skeletonize
    skel = bwskel(imbinarize(smoothed));
end
