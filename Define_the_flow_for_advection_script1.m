
global  U_w U_a U_i V_w V_a V_i eta_x eta_y Th_i T_w A_i%Environmental variables
global  dx dy dt %Grid variables
global  rho_w c_w c_dw rho_a c_a c_da f0 beta g rho_i c_iv T_berg rho_b  %Physical constants.


%Defining the flow of the ocean
for i=1:Nx
    for j=1:Ny
        eta(i,j)=0.008*exp(-((((x(i)-Lx/2)/(Lx/4)).^2))-(((y(j)-(Ly/2))/(Ly/3)).^2));
       % eta(i,j)=0.02*exp(-((((x(i)-Lx/2)/(Lx/4)).^2))-(((y(j)-(Ly/2))/(Ly/3)).^2));
    end
end
for i=1:Nx
    eta_y(i,:)=(eta(i,[2:Ny 1])-eta(i,[Ny 1:(Ny-1)]))./(2*dy);
end
for j=1:Ny
    eta_x(:,j)=(eta([2:Nx 1],j)-eta([Nx 1:(Nx-1)],j))./(2*dx);
end
for i=1:Nx
    for j=1:Ny
        U_w(i,j)=0.0;%(-g/(f0+beta*y(j)).*eta_y(i,j));
        V_w(i,j)=0;%(g/(f0+beta*y(j)).*eta_x(i,j));

    end
end



%Defining the Sea Surface Temperature
for i=1:Nx
    for j=1:Ny
       T_w(i,j)=0;%0.5*exp(-((((x(i)-Lx/2)/(Lx/4)).^2))-(((y(j)-(Ly/2))/(Ly/3)).^2))-1.9;
    end
end
%figure;;pcolor(x/1000,y/1000,T_w');shading('interp');colorbar;


%Defining the flow of the atmosphere
for i=1:Nx
    U_a(i,:)=-1.8;
end
for j=1:Ny
    V_a(:,j)=0.0;%1*cos(4*pi*x/Lx);
end


%Defining the sea ice velcity and thickness, and Area
for i=1:Nx
    for j=1:Ny
        Th_i(i,j)=0;%*exp(-((((x(i)-Lx/5)/(Lx/16)).^2))-(((y(j)-(Ly/5))/(Ly/23)).^2))+4*exp(-((((x(i)-(4*Lx/5))/(Lx/8)).^2))-(((y(j)-((3*Ly)/5))/(Ly/6)).^2));
        U_i(i,j)=0.0;
        V_i(i,j)=0.0;
    end
end
A_i=Th_i./(max(max(Th_i))); % Sea Ice Area

%Defining the land mass 
Land=zeros(Nx,Ny);
Land(1,:)=1;
Land(Nx,:)=1;
Land(:,1)=1;
Land(:,Ny)=1;

Land(Nx,:)=1;


%Land(1:5,:)=1;


x_center=0;
y_center=Lx/2;%3*Ly/6;
x_center2=3*Lx/10;
y_center2=3*Ly/6;


x_center=-Ly/6;2*Ly/3;%5*Ly/14+Lx/9;;
y_center=Ly/2;%3*Ly/6;
x_center2=7*Ly/12;
y_center2=Ly/2;%3*Ly/6;
R_land=Ly/2;
for i=1:Nx
    for j=1:Ny
        
        %Semi circle outwards
%          if (sqrt( (((x(i)-x_center)^2))+(((y(j)-y_center)^2)/0.4)))<R_land
%             Land(i,j)=1;
%          end    
         
          %Semi circle outwards
          if (sqrt( (((x(i)-x_center2)^2))+((y(j)-y_center2)^2)))>R_land   & x(i)<x_center2
            Land(i,j)=1;
         end


    end
end
%Land(2,:)=1;
%Land(Nx-1,:)=1;

x_center2=Lx/2;
y_center2=Ly/2;

eta(:,:)=0;
for i=1:Nx
    for j=1:Ny
        eta(i,j)=0.3*exp(-((((x(i)-x_center2)/(Lx/6)).^2))-(((y(j)-y_center2)/(Ly/6)).^2));
      %  eta(i,j)=0.3*exp(-((((x(i)-Lx/2)/(Lx/8)).^2))-(((y(j)-(Ly/2))/(Ly/8)).^2));
    end
end

for i=1:Nx
    eta_y(i,:)=(eta(i,[2:Ny 1])-eta(i,[Ny 1:(Ny-1)]))./(2*dy);
end
for j=1:Ny
    eta_x(:,j)=(eta([2:Nx 1],j)-eta([Nx 1:(Nx-1)],j))./(2*dx);
end
% for i=1:Nx
%     for j=1:Ny
%         U_a(i,j)=(-g/(f0+beta*y(j)).*eta_y(i,j));
%         V_a(i,j)=(g/(f0+beta*y(j)).*eta_x(i,j));
% 
%     end
% end
%U_w=U_a;
%V_w=V_a;

%Plotting the sea surface height
figure(1);pcolor(x/1000,y/1000,eta');shading('interp');colorbar;
set(gca,'color','k')
caxis([-0 .3]);
xlabel('x axis (km)');
ylabel('y axis (km)');
y_min=0; y_max=Ly/1000; x_min=0;x_max=Lx/1000;
axis([x_min x_max y_min y_max]);drawnow
%pcolor(x/1000,y/1000,Land');shading('interp');colorbar;

hold on
for i=1:Nx
    for j=1:Ny
        if Land(i,j)==1
            plot(x(i)/1000,y(j)/1000,'k*')
            plot((x(i)/1000)+(((dx/2)/1000).*cos(circ_ind)),(y(j)/1000)+(((dx/2)/1000).*sin(circ_ind)),'k','linewidth',4);

        end
    end
end


land_hist=zeros(Nx,Ny);
land_hist(find(Land==1))=1;
