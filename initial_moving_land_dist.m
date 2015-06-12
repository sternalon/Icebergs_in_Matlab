%Initializing Ice_front which is equal to the land
Num_moving_land=sum(sum(Land));  %Matching ice front to all the land.
%Num_moving_land=1;

moving_land_r_n=zeros(Num_moving_land,2);
moving_land_u_n=zeros(Num_moving_land,2);
moving_land_L=zeros(Num_moving_land,1);

% % % %One land particle testing
% moving_land_r_n(1,1)=Lx/3;   moving_land_r_n(1,2)=Ly/2+dx/5;
% moving_land_L(1,1)=dx;% 
% moving_land_u_n(1,1)=0.3;


%Matching moving land to land.
moving_land_L(:,1)=dx;% 
moving_land_u_n(:,1)=0.3;
land_count=0;
for i=1:Nx
    for j=1:Ny
        if Land(i,j)==1
            land_count=land_count+1;
            moving_land_r_n(land_count,1)=x(i);
            moving_land_r_n(land_count,2)=y(j);
            
            if j==1 | j==Ny %| i=Nx
                moving_land_u_n(land_count,:)=0;
            end
        end
    end
end





