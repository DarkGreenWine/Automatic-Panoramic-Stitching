I1=zeros(row,column,3);
I2=I1;
I=I1;
i=1;
theta_i1=phi(4*(i-1)+1);
theta_i2=phi(4*(i-1)+2);
theta_i3=phi(4*(i-1)+3);
R1=expm([0 -theta_i3 theta_i2; theta_i3 0 -theta_i1; -theta_i2 theta_i1 0]);
f1=phi(4*i);
K1=diag([f1,f1,1]);
K1R1=K1*R1;
im1=imread('..\data\building2.jpg');
%reshape(image(:,1),row,column);

i=2;
theta_i1=phi(4*(i-1)+1);
theta_i2=phi(4*(i-1)+2);
theta_i3=phi(4*(i-1)+3);
R2=expm([0 -theta_i3 theta_i2; theta_i3 0 -theta_i1; -theta_i2 theta_i1 0]);
f2=phi(4*i);
K2=diag([f2,f2,1]);
K2R2=K2*R2;
im2=imread('..\data\building3.jpg');
%reshape(image(:,2),row,column);

%theta range 60 to 120 degrees
theta_min=30*pi/180;
theta_max=150*pi/180;
%phi range -90 to 90 degrees
phi1_min=-90*pi/180;
phi1_max=90*pi/180;
%converting to cartitian
%i->theta
%j->phi
for i=1:row
    theta=(theta_max-theta_min)*i/(row)+theta_min;
    for j=1:column
        phi1=(phi1_max-phi1_min)*j/(column)+phi1_min;
        x=sin(theta)*cos(phi1);
        y=sin(theta)*sin(phi1);
        z=cos(theta);
        p=R1*[-z;-y;x];
        %p=[x;y;z];
        xi=floor(p(1)*f1/p(3)+row/2);
        yi=floor(p(2)*f1/p(3)+column/2);
        if((0<xi)&&(xi<row) && (0<yi)&&(yi<column))
            I1(i,j,:)=im1(xi,yi,:);
        end
        p=R2*[-z;-y;x];
        %p=[x;y;z];
        xi=floor(p(1)*f2/p(3)+row/2);
        yi=floor(p(2)*f2/p(3)+column/2);
        if((0<xi)&&(xi<row) && (0<yi)&&(yi<column))
            I2(i,j,:)=im2(xi,yi,:);
        end        
    end
end
W1=(logical(I1)+1)/2;
W2=(logical(I2)+1)/2;
I=(I1+I2)./(W1+W2);
imshow(uint8(I));

