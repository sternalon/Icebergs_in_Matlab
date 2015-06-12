function [accel a_np1 b_np1 a_np2]= accel(berg_count,r_n,u_n,u_star,dt1,Th,l,w,N,M_n,near_a_wall,L_eff,Ranga_not_Velet)
clear u_np1;

global  U_w U_a U_i V_w V_a V_i eta_x eta_y Th_i T_w A_i%Environmental variables
global  dx dy Lx Ly %Grid variables
global  rho_w c_w c_dw rho_a c_a c_da c_wall f0 beta g rho_i c_iv c_ia T_berg rho_b Num_bergs  %Physical constants
global use_periodic_boundaries spring_forces_on relative_damping_on %on off flags
global L_Bond Bond
global F_field
global dt

u0=u_n(berg_count,:);

%Temporary
%u_star=u_n(berg_count,:)+((f0*dt/2)*[u_n(berg_count,2) -u_n(berg_count,1)]); %Crank-Nicholson Coriolis
%u_star=u_n(berg_count,:); %Fully implicit  Coriolis
%u_star=u_star+((f0*dt/2)*[u_n(berg_count,2) -u_n(berg_count,1)]);  %Crank-Nicholson Coriolis

if use_periodic_boundaries==1
    for i=1:Num_bergs
        [r_n(i,1), r_n(i,2)]=Enforce_boundary_conditions(r_n(i,:));
    end
end

exp_coef=2-Ranga_not_Velet;   %exp_coef=1 for RK,  exp_coef=2 for Verlet

%Deciding which velocity to use
r_ind(1)=round((r_n(berg_count,1)/dx))+1;
r_ind(2)=round((r_n(berg_count,2)/dy))+1;

u_w_temp=[U_w(r_ind(1),r_ind(2)) V_w(r_ind(1),r_ind(2))]; %Ocean Velocity.
u_a_temp=[U_a(r_ind(1),r_ind(2)) V_a(r_ind(1),r_ind(2))]; %Atmosphere Velocity
u_i_temp=[U_i(r_ind(1),r_ind(2)) V_i(r_ind(1),r_ind(2))]; % Sea Ice Velocity
Th_i_temp =Th_i(r_ind(1),r_ind(2)); % Sea Ice Thickness

%Finding parameters used in Wave Radiation equation from Martin and Adcroft
a=0.010125*(norm(u_a_temp-u_w_temp).^(2)); %Wave Amplitude
L_w=0.32*(norm(u_a_temp-u_w_temp).^2);  %Wave length
L_t=0.25*L_w; %Upper limit length
L_c=0.125*L_w; %Cuttoff length
c_r=0.06*min(max(0,((l-L_c)/(L_t-L_c))),1); % Wave radiation coeffitient

% Calculate the dimension and mass of iceberg
M=(rho_b*Th*l*w*N);
D=(rho_b/rho_w).*Th; %Draft of iceberg
F=Th-D; %Freeboard of iceberg

%Note that Alon has added a requirement for Bonded Bergs - need to check this.
A_vw1=N*(D-Th_i_temp).*l.*((6-sum(Bond(:,berg_count)))./6);  %Submerged Area of side 1 of the iceberg
A_vw2=N*(D-Th_i_temp).*w.*((6-sum(Bond(:,berg_count)))./6);  %Submerged Area of side 2 of the iceberg
A_va1=N*F.*l.*((6-sum(Bond(:,berg_count)))./6);  %Above water Area of side 1 of iceberg
A_va2=N*F.*w.*((6-sum(Bond(:,berg_count)))./6);  %Above water Area of side 2 of iceberg
A_hw=N*l.*w;  %Area of the base of the iceberg
A_ha=N*l.*w;  %Area of the top of the iceberg


%Calculating force on iceberg - follwing Hunke and Comeau 2011 and Martin and Adcroft 2010

%Non-drag forces are calculated
%F_ss=(M*(f0+beta*r(2))).*[-u_w_temp(2) u_w_temp(1)];%Hunke and Comeau 2011
F_ss=-M*g*[eta_x(r_ind(1),r_ind(2)),eta_y(r_ind(1),r_ind(2))];%Martin and Adcroft 2010
%F_ss=[0 0];
F_r=[0 0];%N*0.5.*rho_w.*c_r.*g.*a*min(a,F).*((2*l*w)./(l+w)).*(u_a_temp./norm(u_a_temp));

%Calculating the force due to interactions with other icebergs
F_ia=F_field(berg_count,:);

c_b=c_wall*near_a_wall; %Ie: high damping is applied when you are near a wall.


Ac=-((f0+beta*r_n(berg_count,2))).*([-u_star(1,2) u_star(1,1)]);

CN_on=1;  %Crank-Nicolson Coriolis on.

p=((f0+beta*r_n(berg_count,2))*((dt/2)));  %Crank-Nicholson Coriolis
%p=((f0+beta*r_n(berg_count,2))*((dt)));     %Fully implicit  Coriolis
if CN_on==0  
    p=p*2;  %Note that Fully implicit  Coriolis does not make inertial circles.
end
a_n_star=(Ac*(CN_on))+((1/(M))*(F_ss+F_r+F_ia)*(~Ranga_not_Velet));  %Things which are half and half.
b_exp=(Ac*(~CN_on))+((1/(M))*(F_ss+F_r+F_ia)*(Ranga_not_Velet));  %Things which are not half half

u_pred=u0;
%Start loop here:
for loop_count=1:2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Relative motion damping
    if relative_damping_on & spring_forces_on==1
        [P_ia P_ia_times_ui]= relative_motion_damping(berg_count,u0,u_pred,r_n,u_n,L_eff,M_n);
    end
    
    %Coefitients adjusted using intermediate step
    C_w=((0.5.*rho_w.*c_w*[A_vw1 A_vw2])+(rho_w.*c_dw*A_hw)).*0.5*(norm(u_w_temp-u0)+norm(u_w_temp-u_pred));
    C_a=((0.5.*rho_a.*c_a*[A_va1 A_va2])+(rho_a.*c_da*A_ha)).*0.5*(norm(u_a_temp-u0)+norm(u_a_temp-u_pred));
    C_i=[0 0];%0.5.*rho_i.*c_iv.*[L W].*Th_i(r_ind(1),r_ind(2)).*0.5*(norm(u_i_temp-u)+norm(u_i_temp-u_pred));
    C_b=(((rho_w.*c_b*A_hw)).*0.5*(norm(-u0)+norm(-u_pred))).*[1 1];
    
    q=(1+((dt/(M)).*(C_w+C_a+C_i+C_b))); %Updated for Verlet
    Factor=[q(1) -p;p q(2)]+((dt/(M))*P_ia);
    Q=Factor^(-1);
    
   % u_np1=(Q*(((u_star')'+((dt/(exp_coef*M))*(F_ss+F_r+F_ia))+((dt/(M))*((C_w.*u_w_temp)+(C_a.*u_a_temp)+(C_i.*u_i_temp)+C_ia_times_ui+P_ia_times_ui)))'))';    
   % AX=(1/dt)*(Q*((dt*Ac+((dt/(exp_coef*M))*(F_ss+F_r+F_ia))+((dt/(M))*((C_w.*(u_w_temp-u_star))+(C_a.*(u_a_temp-u_star))+(C_i.*(u_i_temp-u_star))+(P_ia_times_ui-((P_ia*u_star')')))))'))';
   % AX=(Q*((0.5*Ac+((1/(exp_coef*M))*(F_ss+F_r+F_ia))+((1/(M))*((C_w.*(u_w_temp-u_star))+(C_a.*(u_a_temp-u_star))+(C_i.*(u_i_temp-u_star))+(P_ia_times_ui-((P_ia*u_star')')))))'))';
     AX=(Q*(((a_n_star/2)+b_exp+((1/(M))*((C_w.*(u_w_temp-u_star))+(C_a.*(u_a_temp-u_star))+(C_i.*(u_i_temp-u_star))+(P_ia_times_ui-((P_ia*u_star')')))))'))';
    u_pred=u_star+(dt*AX);    
end
u_np1=u_pred;

% %Now, the full forces can be defined
% F_w=C_w.*(u_w_temp-u_np1);
% F_a=C_a.*(u_a_temp-u_np1);
% F_i=[0 0];%C_i.*(u_i_temp-u_np1);
% F_c=-(M*0.5*(f0+beta*r_n(berg_count,2))).*([-u_np1(1,2) u_np1(1,1)]);  %Crank-Nicholson
% %F_c=-(M*(f0+beta*r_n(berg_count,2))).*([-u_np1(1,2) u_np1(1,1)]);  %Fully implicit coriolis
% %F_c=-(M*0.5*(f0+beta*r_n(berg_count,2))).*([-u_np1(1,2) u_np1(1,1)]+[u_n(berg_count,2) -u_n(berg_count,1)]);  %Crank-Nicholson added inside
% F_b=C_b.*(-u_np1);
% F_ia_damp=P_ia_times_ui-((P_ia*(u_np1'))');


Ac=-((f0+beta*r_n(berg_count,2))).*([-u_np1(1,2) u_np1(1,1)]);
a_np1=(Ac*(CN_on))+((1/(M))*(F_ss+F_r+F_ia)*(~Ranga_not_Velet));  %Things which are half and half.
accel=AX;
b_np1=AX-(a_np1/2);   %?????

% 
% accel=(1./M).*(F_ss+F_r+F_ia+F_w+F_a+F_c+F_i+F_ia_damp+F_b);
% a_np1=(1./M).*(F_ss+F_r+F_ia);
% b_np1=(1./M).*(F_w+F_a+F_c+F_i+F_ia_damp+F_b);



