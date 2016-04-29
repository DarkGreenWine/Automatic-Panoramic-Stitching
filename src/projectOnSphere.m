% This function will take an image 
% And the registration parameters
% Convert it to spherical system
function [ out ] = projectOnSphere(I,f)
%I = imread('building1.JPG');
%f = 400;
width = size(I,2);
height = size(I,1);

% Calculate image center. Middle pixel.
x_center = floor(width/2);
y_center = floor(height/2);

unwrapped = zeros(floor(f), floor(f), 3);

for i = 1:width
    for j = 1:height
        curr_x = i-x_center; 
        curr_y = j-y_center; 

	    mag_XYZ = sqrt(curr_x^2 + curr_y^2 + f^2);
        phi = asin(curr_y/mag_XYZ);
        theta = asin(curr_x/(mag_XYZ*cos(phi)));

        x_sph = round(f*theta + x_center);
        y_sph = round(f*phi + y_center);
        unwrapped(y_sph,x_sph,:) = I(j,i,:);
    end
end
out = uint8(unwrapped);
imshow(out);















