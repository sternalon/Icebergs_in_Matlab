clear all
close all

load 'blue_white.mat'

path_name='/Users/alon/Desktop/files/Icebergs_clusters/Matlab_scripts/Box_experiment';
cd(path_name)



%filename='Viscous_spring_contact'
%filename='Different_sizes_Viscous_spring_contact'
%filename='Including_big_Bergs_with_Jacobi_5';
%filename='All_bonds_Including_big_Bergs_with_Jacobi_5';
%filename='Half_circle_exp_with_bergs_Jacobi_5';
%filename='Half_circle_exp_with_little_bergs_Jacobi_5';
%filename='Two_interacting_bergs_springs_for_land_bonds_interaction';
filename='Test_using_L_eff';
%filename='Test_using_L_eff_alpha_max_10';
%filename='Box_using_L_eff_alpha_max_10';
%filename='Conc_on_Box_using_L_eff_alpha_max_10';
%filename='Convex__Conc_on_alpha_max_10';
%filename='Convex__Conc_on_alpha_mu_0';
%filename='Convex__Conc_on_alpha_mu_1';
%filename='Convex__Conc_on_alpha_mu_3';
%filename='Convex_interaction_factor_mu_3';
%filename='Convex_interaction_factor_mu_3_gamma_0_3';
%filename='Convex_interaction_factor_mu_3_gamma_0_3_take2';
%filename='Convex_interaction_factor_mu_3_gamma_0_3_take3';
%filename='Convex_interaction_factor_mu_5_gamma_1_take4';
filename='Convex_interaction_factor_mu_5_gamma_1_take5';
%filename='Reverse_wind_dir';
%filename='Circular_wind_dir';
%filename='Sea_surface_slope';
%filename='Gridded_bergs_towards_wall';
% %filename='Corriolis_Gridded_bergs_towards_wall';
% %filename='Gridded_bergs_towards_wall_again';
% filename='Gridded_bergs_towards_wall_again2';
% filename='Gridded_bergs_towards_wall_again_stiffer';
% filename='Gridded_bergs_towards_wall_again_stiffer_2';
% %filename='Gridded_bergs_towards_wall_again_stiffer_Jacobi';
% filename='Gridded_bergs_towards_wall_again_stiffer_Jacobi';
% filename='Gridded_bergs_towards_wall_again_stiffer_Jacobi_verlet';
% filename='Gridded_bergs_towards_wall_again_stiffer_Jacobi_verlet_with_contact_correction';
 filename='Two_particle_spring';
% filename='Two_particle_Jacobi';
 filename='Big_and_small_bergs_curved_coast_again';
%  filename='New_test_curved_coast';
% % filename='New_test_curved_coast_up_pia';
%filename='New_test_curved_coast_up_pia_and_q_ia';
% filename='Test_new_forces';
% filename='Footloose_test1_Linit_100km'
% filename='Footloose_test1_Linit_100km_movie';
% filename='Footloose_test1_Linit_100km_Movie2';
% %filename='From_saved_test';
% filename='Calving_icebergs1';
% filename='Calving_icebergs2';
% filename='Calving_icebergs3';
% filename='Calving_icebergs3';
 filename='Calving_icebergs4';
 filename='Calving_icebergs5';
 filename='Calving_icebergs6';

circ_ind=[0:0.1:2*pi];
init_time=40;
time_step=1;
final_time_ind=0;   %final_time=0 to automatically do until end of list

plot_the_area=0;
plot_dots_and_circ=1;
plot_the_bonds=0; 
plot_land_each_time=0;
plot_L_eff=1;
make_a_movie=1;

%Defining the grid
Lx=500*1000;%100*1000; %Length of the domain
Ly=500*1000;%120*1000 %Width of the domain
dx=10*1000;  %ocean grid spacing in x
dy=dx; %ocean grid spacing in y
Nx=(Lx/dx)+1; %number of positions in the x grid
Ny=(Ly/dy)+1; %number of positions in the y grid
x=([1:Nx]-1).*dx;
y=([1:Ny]-1).*dy;
fontsize=18;


Hov_Conc=zeros(302,101);

cd data
cd(filename)

load 'iceberg_mat_data.mat'  %berg_pos_x berg_pos_y berg_pos_L berg_pos_time
load 'environmental_variables.mat'  %Land eta U_w V_w x y Num_bergs
load 'iceberg_time_ind.mat'
load 'Parameters.mat'


%Defining the grid for concentration
dx=2*1000;  %ocean grid spacing in x
dy=dx;
Nx2=(Lx/dx)+1; %number of positions in the x grid
Ny2=(Ly/dy)+1; %number of positions in the y grid
x2=([1:Nx2]-1).*dx;
y2=([1:Ny2]-1).*dy;
W=berg_pos_L(:,1);

if final_time_ind==0;
    final_time_ind=length(berg_pos_time);
end

time_ind_list=[init_time:time_step:final_time_ind];


if make_a_movie==1
    'Begining to write the movie'
    addpath(genpath('/Users/alon/Desktop/files/Random/export_fig/'));
    nframes = length(time_ind_list);
    framerate =8;% 10;    % in frames per second, for movie making.
    writerObj = VideoWriter(strcat(filename,'.avi'));
    writerObj.FrameRate = framerate;
    open(writerObj);

end
cd ../
cd ../


%figure;
figure('units','normalized','outerposition',[0 0 1 1]);
%axes1=axes('Position',[0.2 0.2 0.99 0.99])
%Plotting the initial setup
if plot_dots_and_circ==1
    time_count=time_ind_list(1);
    time=berg_pos_time(time_count);
    % hold on
    for berg_count=1:Num_bergs
        hline(:,berg_count)=plot(berg_pos_x(berg_count,time_count)/1000,berg_pos_y(berg_count,time_count)/1000,'*r');
        hline2(:,berg_count)=plot((berg_pos_x(berg_count,time_count)/1000)+(((berg_pos_L(berg_count,time_count)/2)/1000).*cos(circ_ind)),(berg_pos_y(berg_count,time_count)/1000)+(((berg_pos_L(berg_count,time_count)/2)/1000).*sin(circ_ind)),'g');
        if plot_L_eff==1
            hline8(:,berg_count)=plot((berg_pos_x(berg_count,time_count)/1000)+(((berg_pos_L_eff(berg_count,time_count)/2)/1000).*cos(circ_ind)),(berg_pos_y(berg_count,time_count)/1000)+(((berg_pos_L_eff(berg_count,time_count)/2)/1000).*sin(circ_ind)),'r');
        end
        hline3=text((Lx/2/1000)-(Lx/20/1000),(Ly/1000)+(Ly/20/1000),['time = ' num2str(round(10*time/60/60/24)/10) ' days'],'fontsize',fontsize);
    end
end



pcolor(x/1000,y/1000,eta');shading('interp');%colorbar;
set(gca,'color','k')
caxis([-0 1]);
xlabel('x axis (km)');
ylabel('y axis (km)');




%Plotting the land
count=0;
hold on
for i=1:Nx
    for j=1:Ny
        if Land(i,j)==1
            count=count+1;
            hline10(i,j)=plot(x(i)/1000,y(j)/1000,'k*');
            hline11(i,j)=plot((x(i)/1000)+(((dx/2)/1000).*cos(circ_ind)),(y(j)/1000)+(((dx/2)/1000).*sin(circ_ind)),'k','linewidth',4);
        end
    end
end


%axis([0 Lx/1000 0 Ly/1000]);drawnow
%max_lat=Ly/1000;min_lat=0;max_lon=Lx/1000;min_lon=0;
%max_lat=Ly/1000;min_lat=0;max_lon=500;min_lon=0;  %Big and small bergs

%max_lat=385;min_lat=15;max_lon=850;min_lon=0;  %Calving 4
max_lat=385;min_lat=15;max_lon=850;min_lon=50;  %Calving 5

% max_lat=400; max_lon=500; min_lat=150; min_lon=10;  %footloose

% max_lat=350; max_lon=700; min_lat=150; min_lon=400;  %Two particle spring
 axis([min_lon max_lon min_lat max_lat]);drawnow

%axis([400 900 0 Ly/1000]);drawnow
%axis([0 500 0 500]);drawnow
%axis([450 650 150 350]);drawnow

%axis([75 500 0 Ly/1000]);drawnow
hold on
colormap(blue_white)
for time_count=time_ind_list;
    time=berg_pos_time(time_count);
    
    if plot_the_area==1
        Conc=zeros(Nx2,Ny2);
        %    Calculating area on a grid
        if time_count>init_time %& berg_count==1
            delete(hline5);
            if plot_land_each_time==1
                delete(hline10);
                delete(hline11);
            end
            
            % delete
        end
        for i=1:Nx2
            for j=1:Ny2
                alpha=3;
                R=sqrt(  ((berg_pos_x(:,time_count)-x2(i)).^2) +  ((berg_pos_y(:,time_count)-y2(j)).^2));
                Conc(i,j)=sum(pi.*((berg_pos_L(:,time_count)./2).^(2)).*Kernal(R,berg_pos_L_eff(:,time_count),alpha));
                
                % R=sqrt(  ((berg_pos_x(:,time_count)-x(i)).^2) +  ((berg_pos_y(:,time_count)-y(j)).^2));
                % Conc(i,j)=sum(pi.*(L(:,1)./2).^(2).*Kernal(R,L,alpha));
            end
        end
        
        
        if max(max(Conc))>1.1
            ['There is a very large Conc at time = ' num2str(round(10*time/60/60/24)/10) ' days']
        end
        %Conc(find(Land==1))=NaN;
        hline5=pcolor(x2/1000,y2/1000,Conc');shading('interp');caxis([0. 1]);%colorbar
        %colormap gray
        
        %Plotting the land
        if plot_land_each_time==1
            hold on
            for i=1:Nx
                for j=1:Ny
                    if Land(i,j)==1
                        hline10(i,j)=plot(x(i)/1000,y(j)/1000,'k*');
                        hline11(i,j)=plot((x(i)/1000)+(((dx/2)/1000).*cos(circ_ind)),(y(j)/1000)+(((dx/2)/1000).*sin(circ_ind)),'k','linewidth',4);
                        
                    end
                end
            end
        end
    end
    
    if plot_dots_and_circ==1
        for berg_count=1:Num_bergs
            if time_count>0 %& berg_count==1
                delete(hline(:,berg_count));
                delete(hline2(:,berg_count));
                if plot_L_eff==1
                    delete(hline8(:,berg_count));
                end
            end
           hline(:,berg_count)=plot(berg_pos_x(berg_count,time_count)/1000,berg_pos_y(berg_count,time_count)/1000,'.w');
         %   hline(:,berg_count)=plot(berg_pos_x(berg_count,time_count)/1000,berg_pos_y(berg_count,time_count)/1000,'.w','MarkerSize',260);%two particle spring
         
        %    hline(:,berg_count)=plot(berg_pos_x(berg_count,time_count)/1000,berg_pos_y(berg_count,time_count)/1000,'.w','MarkerSize',90);  
           % hline2(:,berg_count)=plot((berg_pos_x(berg_count,time_count)/1000)+(((berg_pos_L(berg_count,time_count)/2)/1000).*cos(circ_ind)),(berg_pos_y(berg_count,time_count)/1000)+(((berg_pos_L(berg_count,time_count)/2)/1000).*sin(circ_ind)),'g');
            hline2(:,berg_count)=fill((berg_pos_x(berg_count,time_count)/1000)+(((berg_pos_L(berg_count,time_count)/2)/1000).*cos(circ_ind)),(berg_pos_y(berg_count,time_count)/1000)+(((berg_pos_L(berg_count,time_count)/2)/1000).*sin(circ_ind)),'w');

            if plot_L_eff==1
                hline8(:,berg_count)=plot((berg_pos_x(berg_count,time_count)/1000)+(((berg_pos_L_eff(berg_count,time_count)/2)/1000).*cos(circ_ind)),(berg_pos_y(berg_count,time_count)/1000)+(((berg_pos_L_eff(berg_count,time_count)/2)/1000).*sin(circ_ind)),'w');
            end
        end
    end
    
    if plot_the_bonds==1
        %Plotting the bonds
        for berg_count=1:Num_bergs
            for pair_berg=1:berg_count
                if (Bond(pair_berg,berg_count)==1)
                    if time_count>init_time
                        delete(hline4(:,berg_count,pair_berg));
                    end
                    hline4(:,berg_count,pair_berg)=plot([berg_pos_x(berg_count,time_count)/1000 berg_pos_x(pair_berg,time_count)/1000],[berg_pos_y(berg_count,time_count)/1000 berg_pos_y(pair_berg,time_count)/1000],'m','linewidth',3);
                end
            end
            if (L_Bond(pair_berg,berg_count)==1)
                if time_count>1
                    delete(hline5(:,berg_count,pair_berg));
                end
                hline5(:,berg_count,pair_berg)=plot([berg_pos_x(berg_count,time_count)/1000 berg_pos_x(pair_berg,time_count)/1000],[berg_pos_y(berg_count,time_count)/1000 berg_pos_y(pair_berg,time_count)/1000],'r','linewidth',3);
            end
        end
    end
    
    
    
    if time_count>0 %& berg_count==1
        delete(hline3);
    end
    %hline3=text((Lx/2/1000)-(Lx/20/1000),(Ly/1000)+(Ly/20/1000),['time = ' num2str(round(10*time/60/60/24)/10) ' days'],'fontsize',fontsize);
    hline3=text((((max_lon-min_lon)/2))-(((max_lon-min_lon)/2)/10),(max_lat)+(max_lat/40),['time = ' num2str(round(10*time/60/60/24)/10) ' days'],'fontsize',fontsize);
    set(gca,'fontsize',fontsize) 
    drawnow
    
    if make_a_movie==1
        %Using regular method:
        %frame = getframe;
        
        %Using Mitch's method with non-visible figures
        set(gcf,'visible','off') 
        frame = export_fig(gcf,'-nocrop', '-a1');
        
        writeVideo(writerObj,frame);
        ['Writing movie frame at time t=' num2str(time/(60*60*24)) ' days']
    end
    
end

if make_a_movie==1
    %cd data
    %cd(filename)
    close(writerObj);
    'Writing movie done!'
    %cd ..
    %cd ..
end





