%Saves still frames from start and end of .tiff stack, to be processed for
%masking

set(0,'DefaultFigureWindowStyle','docked')


%load image start
I_start = imread(data,'Index',start);
t = imgaussfilt(I_start(:,:,1),[50 50]);
t2 = (I_start(:,:,1)-t) == 0; 
t3 = medfilt2(t2,[20 20]);
t4_start = ~t3;
figure(); imagesc(t4_start);

%load image end
I_end = imread(data,'Index',numFrames);
t = imgaussfilt(I_start(:,:,1),[50 50]);
t2 = (I_end(:,:,1)-t) == 0; 
t3 = medfilt2(t2,[20 20]);
t4_end = ~t3;
figure(); imagesc(t4_start);

%save images
imwrite(t4_start,[source '_mask_ref_start.png']);
imwrite(t4_end,[source '_mask_ref_end.png']);