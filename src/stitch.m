%function [out] = stitch(img1,img2,f)
f = 100;
img1 = imread('../data/building1.JPG');
img2 = imread('../data/building2.JPG');
sph_img1 = projectOnSphere(img1,f);
sph_img2 = projectOnSphere(img2,f);

% The sph_image contains black parts
x = size(sph_img1,1);
y = size(sph_img1,2);

out1 = remove_black(sph_img1);
out2 = remove_black(sph_img2);
out = [out1,out2];

