n=2;
im=imread('..\data\building2.jpg');
im=im(:,:,1);
[row,column]=size(im);
image=zeros(row*column,n);
image(:,1)=im(:);

im=imread('..\data\building3.jpg');
im=im(:,:,1);
image(:,2)=im(:);

% im=imread('..\data\building3.jpg');
% im=im(:,:,1);
% image(:,3)=im(:);
% 
% im=imread('..\data\building4.jpg');
% im=im(:,:,1);
% image(:,4)=im(:);
% 
% im=imread('..\data\building5.jpg');
% im=im(:,:,1);
% image(:,5)=im(:);
%calculating nearest matches
matching_matrix=zeros(n);
for i =1:n-1
    [ des1, loc1] = sift(reshape(image(:,i),row,column));
    for j=i+1:n
        [ des2, loc2] = sift(reshape(image(:,j),row,column));
        [p1,p2] = match(des1, loc1,des2, loc2);
        matching_matrix(i,j)=length(p1);
        matching_matrix(j,i)=length(p1);
    end
end
% mapping is of size n x nearest_neighbour ith row contains indices of k
% nearest_neighbour images
nearest_neighbour=1;
[~,mapping]=sort(matching_matrix,2,'descend');
mapping=mapping(:,[1:nearest_neighbour]);

K=4; %for homography
mapping_points=zeros(n*K,nearest_neighbour*4);
for i =1:n
    [ des1, loc1] = sift(reshape(image(:,i),row,column));
    for j=1:nearest_neighbour
        [ des2, loc2] = sift(reshape(image(:,mapping(i,j)),row,column));
        [p1,p2] = match(des1, loc1,des2, loc2);
        [~, inliers] = ransacfithomography(p1', p2', 0.005);
        mapping_points((i-1)*K+1:i*K,(j-1)*4+1:j*4)=[ p1(inliers(1:K),[1,2]),  p2(inliers(1:K),[1,2]) ];
    end
end

