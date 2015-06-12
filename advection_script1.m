%This is a first test for making a particle advection
clear all
close all

global  U_w U_a U_i V_w V_a V_i eta_x eta_y Th_i T_w A_i%Environmental variables.
global  dx dy dt Lx Ly Nx Ny%Grid variables
global  rho_w c_w c_dw rho_a c_a c_da f0 beta g rho_i c_iv T_berg rho_b  c_ia c_wall p_ia %Physical constants
global r_n Num_bergs Bond L_Bond W_max %These are being made global for now, but might change later.
global use_periodic_boundaries Land x y
global Ks C_damp Ks_land C_damp_land Ks_bond C_damp_bond B1 spring_forces_on relative_damping_on use_land_force
global PE_n
global F_field

%% On off Flags
Ranga_not_Velet=0;  %Range=1, Verlet=0;

plot_each_time=1;
write_iceberg_data=0; %For writing the data.
plot_the_area=0;
thermodynamics_on=0;
plot_the_bonds=0;
plot_land_each_time=0;
plot_land_history=0;
plot_movgin_land=0;
plot_just_velocities=0; berg_view_num=1;  %Plots the velocity of one of the bergs.
plot_spring_energy=0; berg_view_num=1; plot_the_total_energy=1;

load_from_file=0;

Jacobi_correction=1;
Jacobi_initial_bond_correction=1;
Jacobi_correct_land=1;
Jacobi_bond_correction=0;
Jacobi_contact_correction=0;
use_moving_land=0;

use_land_force=0;
spring_forces_on=1;
relative_damping_on=1;

Rolling_icebergs=0;
Turn_all_damping_off=0;
use_periodic_boundaries=0;

one_particle_testing=0;
two_particle_testing=0;
three_particle_testing=0;
three_particle_testing2=0;
three_particle_testing_bonded=0;

%% Initializing variables.

%Calling script to initialize parameters.
parameters_advection_script1


%Defining the ocean grid
Nx=(Lx/dx)+1; %number of positions in the x grid
Ny=(Ly/dy)+1; %number of positions in the y grid
x=([1:Nx]-1).*dx;
y=([1:Ny]-1).*dy;
X=diag(x)*ones(Nx,Ny);
Y=ones(Nx,Ny)*diag(y);

%Initializing the envioronmental variables
U_w=zeros(Nx,Ny);
V_w=zeros(Nx,Ny);
U_a=zeros(Nx,Ny);
V_a=zeros(Nx,Ny);
U_i=zeros(Nx,Ny);
V_i=zeros(Nx,Ny);
eta=zeros(Nx,Ny);
eta_x=zeros(Nx,Ny);
eta_y=zeros(Nx,Ny);
Th_i=zeros(Nx,Ny);
Conc=zeros(Nx,Ny);

%using for visualizatoin
circ_ind=[0:0.1:2*pi];

%Script to Defining the flow of the ocean, atmosphere...
Define_the_flow_for_advection_script1

%Initializing iceberg position
if use_moving_land==1
    initial_moving_land_dist
end

%Defining the initial iceberg distribution
initial_berg_distribution

if load_from_file==1
    cd data
    cd(load_filename)
    load 'iceberg_mat_data.mat'
    load 'environmental_variables.mat'
    load 'Parameters.mat'
    [temp load_time_ind]=size(berg_pos_x);
    load_time_ind=1201;
    
    %Time dependent variables
    r_n(:,1)=berg_pos_x(:,load_time_ind);
    r_n(:,2)=berg_pos_y(:,load_time_ind);
    u_n(:,1)=berg_vel_x(:,load_time_ind);
    u_n(:,2)=berg_vel_y(:,load_time_ind);
    L(:,1)=berg_pos_L(:,load_time_ind);
    L_eff(:,1)=berg_pos_L_eff(:,load_time_ind);
    % l(:,1)=berg_pos_l(:,load_time_ind);
    % w(:,1)=berg_pos_w(:,load_time_ind);
    % N(:,1)=berg_pos_N(:,load_time_ind);
    % Th(:,1)=berg_pos_Th(:,load_time_ind);
    
    cd ../
    cd ../
    parameters_advection_script1
    
end

if spring_forces_on==0
    c_ia=0.00;
    p_ia=0.00;
end

Testing_senarios

if Turn_all_damping_off==1
    c_w=0;c_dw=0;
    c_a=0;c_da=0;
end

%% More initializing

%Finding initial iceberg index
r_ind(:,1)=round((r_n(:,1)/dx))+1;
r_ind(:,2)=round((r_n(:,2)/dy))+1;

% %Giving the icebergs initial velocity equal to the ocean, when loaded
% if load_from_file==0
%     for berg_count=1:Num_bergs
%         u_n(berg_count,1)=0.0;
%         u_n(berg_count,2)=0.0;
%     end
% end

%Temporarily here for initializing time step
%u_np1=u_n;

%
% %Making the first length bigger than the second
% for berg_count=1:Num_bergs
%     if L(berg_count,1)<W(berg_count,1)
%         temp=L(berg_count,1);
%         L(berg_count,1)=W(berg_count,1);
%         W(berg_count,1)=temp;
%     end
% end

% Calculate the dimension and mass of iceberg
M(:,1)=rho_b*L(:,1).*W(:,1).*Th(:,1); %Mass of iceberg
%D(:,1)=(rho_b/rho_w).*Th(:,1); %Draft of iceberg
%F(:,1)=Th(:,1)-D(:,1); %Freeboard of iceberg

%Defining the initial effective width.
%L_eff=1.5.*L;
%L_eff(find(sum(Bond)>0))=L(find(sum(Bond)>0));

%L_eff(find(N==1))=L(find(N==1));

%% Initializing figures

%figure(1)
hold on
for berg_count=1:Num_bergs
    hline(:,berg_count)=plot(r_n(berg_count,1)/1000,r_n(berg_count,2)/1000,'*r');
    hline2(:,berg_count)=plot((r_n(berg_count,1)/1000)+(((L(berg_count,1)/2)/1000).*cos(circ_ind)),(r_n(berg_count,2)/1000)+(((L(berg_count,1)/2)/1000).*sin(circ_ind)),'g');
    hline8(:,berg_count)=plot((r_n(berg_count,1)/1000)+(((L_eff(berg_count,1)/2)/1000).*cos(circ_ind)),(r_n(berg_count,2)/1000)+(((L_eff(berg_count,1)/2)/1000).*sin(circ_ind)),'m');
end
if use_moving_land==1
    for land_count=1:Num_moving_land
        hline70(:,land_count)=plot(moving_land_r_n(land_count,1)/1000,moving_land_r_n(land_count,2)/1000,'*m');
        hline71(:,land_count)=plot((moving_land_r_n(land_count,1)/1000)+(((moving_land_L(land_count,1)/2)/1000).*cos(circ_ind)),(moving_land_r_n(land_count,2)/1000)+(((moving_land_L(land_count,1)/2)/1000).*sin(circ_ind)),'g');
    end
end

%y_min=160; y_max=340; x_min=20;x_max=160;
%y_min=(Ly/2-(2*dx))/1000; y_max=(Ly/2+(2*dx))/1000; x_min=(Lx/2-(2*dx))/1000;x_max=(Lx/2+(2*dx))/1000;
y_min=0; y_max=Ly/1000; x_min=0;x_max=Lx/1000;
axis([x_min x_max y_min y_max]);drawnow
%axis([0 Lx/1000 0 Ly/1000]);drawnow
grid on

%% Start of the time loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Start of time loop.   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_count=0;
check_total=0;
for time=start_time:dt:end_time;
    time_count=time_count+1;
    
    %Check if particles are near a wall, and if so, increase damping
    near_a_wall=zeros(Num_bergs,1);
    for count=1:Num_bergs
        L_crit=2*(L_eff(count,1)+dx)./2;
        Num_points=ceil(L_crit/dx);
        r_ind=round(r_n(count,:)/dx)+1;
        for q=r_ind(1)-Num_points:r_ind(1)+Num_points
            if (q<(Nx+1)) & (q>0)
                for j=r_ind(2)-Num_points:r_ind(2)+Num_points
                    if (j<(Ny+1)) & (j>0)
                        if Land(q,j)==1
                            %Find out if there is any overlap
                            if (abs((r_n(count,1)-x(q)))<L_crit)  & (abs((r_n(count,2)-y(j)))<L_crit)
                                near_a_wall(count)=1;
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    if Ranga_not_Velet==1
        Ranga_Katta_stepping
    end
    if Ranga_not_Velet==0
        Verlet_stepping;
    end
    
    
    %Manually Correct the positions of icebergs over land and overlapping
    if Jacobi_correction==1
        for k=1:1
            jacobi_iteration;
        end
    end
    
    
    %Enforcing periodic boundary conditions
    if use_periodic_boundaries==1
        for berg_count=1:Num_bergs
            [r_np1(berg_count,1), r_np1(berg_count,2)]=Enforce_boundary_conditions(r_np1(berg_count,:));
        end
    end
    
    %Switching old and new variables
    r_n=r_np1;
    u_n=u_np1;
    
    %Moving the moving land
    if use_moving_land==1
        for land_count=1:Num_moving_land
            moving_land_r_n(land_count,1)=moving_land_r_n(land_count,1)+dt*moving_land_u_n(land_count,1);
            moving_land_r_n(land_count,2)=moving_land_r_n(land_count,2)+dt*moving_land_u_n(land_count,2);
        end
    end
    
       
    %Thermodynamics for melting of the icebergs
    if thermodynamics_on==1
        apply_thermodynamics
    end
    
    %Plots various things.
    run_plotting_script;
    
    if write_iceberg_data==1
       write_data
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF TIME STEPPING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%
%% Write to disk.
if write_iceberg_data==1
    cd data
    cd(filename)
    save 'iceberg_mat_data.mat' berg_pos_x berg_pos_y berg_vel_x berg_vel_y berg_pos_L berg_pos_time berg_pos_L_eff
    save 'iceberg_time_ind.mat' Bond L_Bond Num_bergs
    save 'environmental_variables.mat' Land eta U_w V_w x y dx dy Nx Ny Lx Ly
    save 'Parameters.mat' Omega R_earth g lat_ref f0 beta W_max Th_max T_berg Max_Num Ks C_damp Ks_land C_damp_land Ks_bond C_damp_bond alpha_max C1 gamma mu B1 c_w c_dw rho_a c_a c_da rho_i c_iv c_wall c_ia p_ia Ks_A
    
    cd ../
    cd ../
end





