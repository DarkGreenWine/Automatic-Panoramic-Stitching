z1=[0; 0 ; 0 ; 1]*ones(1,n);
phi=[0;0;0.1;300;0;0;-0.1;300];
phi_prev=zeros(4*n,1);
m=0; 
% v=1.0066;
% lambda=v*271;
% v=1.002;%/1.0009;
% lambda=v*831;
v=1/1.1;%/1.0009;
lambda=v*0.11;
while m<35;%320%sum(abs(phi-phi_prev).^2)>10 % need to change this value
    %calculating JtJ and Jtr matrix
    JtJ=zeros(4*n);
    Jtr=zeros(4*n,1);
    e=0;
    for i=1:n
        % calculating Ri and Ki
        theta_i1=phi(4*(i-1)+1);
        theta_i2=phi(4*(i-1)+2);
        theta_i3=phi(4*(i-1)+3);
        Ri=expm([0 -theta_i3 theta_i2; theta_i3 0 -theta_i1; -theta_i2 theta_i1 0]);
        fi=phi(4*i);
        Ki=diag([fi,fi,1]);
        KiRi=Ki*Ri;
        
        for l=1:nearest_neighbour
            j=mapping(i,l);
            % calculating inverse of Rj and Kj
            theta_j1=phi(4*(j-1)+1);
            theta_j2=phi(4*(j-1)+2);
            theta_j3=phi(4*(j-1)+3);
            RjT=expm([0 theta_j3 -theta_j2; -theta_j3 0 theta_j1; theta_j2 -theta_j1 0]);
            fj=phi(4*j);
            KjInv=diag([1/fj,1/fj,1]);
            RjT_KjInv=RjT*KjInv;
            
            %calculating  the derivatives of H(transformation) wrt phi_i
            H=KiRi*RjT_KjInv;
            Hi1=KiRi*[0 0 0; 0 0 -1; 0 1 0]*RjT_KjInv;
            Hi2=KiRi*[0 0 1; 0 0 0 ;-1 0 0]*RjT_KjInv;
            Hi3=KiRi*[0 -1 0; 1 0 0; 0 0 0]*RjT_KjInv;
            Hif=diag([1 1 0])*Ri*RjT_KjInv;
            
            %calculating  the derivatives of H(transformation) wrt phi_j
            KiRiRjT=KiRi*RjT;
            Hj1=KiRiRjT*[0 0 0; 0 0 1; 0 -1 0]*KjInv;
            Hj2=KiRiRjT*[0 0 -1; 0 0 0 ;1 0 0]*KjInv;
            Hj3=KiRiRjT*[0 1 0; -1 0 0; 0 0 0]*KjInv;
            Hjf=KiRiRjT*diag([-1/fj^2 -1/fj^2 0]);
            %Hjfextra=H*[1; 0; 0];
            %Hjfextra=H*[0; 0; 1];
            for k=1:4
                uj=[mapping_points(4*(i-1)+ k ,4*(l-1)+3:4*l) , 1]'; 
                %uj=[fj;-uj(2);-uj(1)];
                ui=mapping_points(4*(i-1)+ k ,4*(l-1)+1:4*(l-1)+2)';
                pij_k=H*uj;
                %rij_k=ui-[-pij_k(3);-pij_k(2)];
                rij_k=ui-pij_k(1:2);
                x=pij_k(1);
                y=pij_k(2);
                z=pij_k(3);
                
                dri=[1/z 0 -x/z^2 ;0 1/z -y/z^2]*[Hi1*uj Hi2*uj Hi3*uj Hif*uj];
                drj=[1/z 0 -x/z^2 ;0 1/z -y/z^2]*[Hj1*uj Hj2*uj Hj3*uj Hjf*uj];
                %ii
                JtJ(4*(i-1)+1:4*i, 4*(i-1)+1:4*i)=JtJ(4*(i-1)+1:4*i, 4*(i-1)+1:4*i)+dri'*dri;
                %ij               
                JtJ(4*(i-1)+1:4*i, 4*(j-1)+1:4*j)=JtJ(4*(i-1)+1:4*i, 4*(j-1)+1:4*j)+dri'*drj;
                %ji
                JtJ(4*(j-1)+1:4*j, 4*(i-1)+1:4*i)=JtJ(4*(j-1)+1:4*j, 4*(i-1)+1:4*i)+drj'*dri;
                %jj
                JtJ(4*(j-1)+1:4*j, 4*(j-1)+1:4*j)=JtJ(4*(j-1)+1:4*j, 4*(j-1)+1:4*j)+drj'*drj;
                
                %ith
                Jtr(4*(i-1)+1:4*i)=Jtr(4*(i-1)+1:4*i)+dri'*rij_k;
                %jth
                Jtr(4*(j-1)+1:4*j)=Jtr(4*(j-1)+1:4*j)+drj'*rij_k;
                e=e+sum(rij_k.^2);
            end
        end
    end
    
    sigma_t=pi/16;
    favg=sum(phi(4*(1:n)))/n;
    sigma_f=favg/10;
    z1=[sigma_t^2; sigma_t^2 ; sigma_t^2 ; sigma_f^2]*ones(1,n);
    Cp=diag(z1(:));
    lambda=lambda/v;
    phi_prev=phi;
    phi=phi+(JtJ+lambda*diag(diag(JtJ)))\Jtr;
    m=m+1;
    
end
e