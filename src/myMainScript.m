%% MyMainScript

% Read all the images
% the folder in which the images exists
source = uigetdir();
source1 = strcat(source,'/*.jpg');
srcImages = dir(source1);  
% img is (x,y,N) matrix
for i = 1 : length(srcImages)
    filename = strcat(source,'/',srcImages(i).name);
    img(:,:,i) = imread(filename);
    %figure, imshow(img(:,:,i);
end

% Extract SIFT features from all n images

