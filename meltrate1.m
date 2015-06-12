function [M_b , M_e ,M_v]= meltrate1(r,u,L)

global  U_w U_a V_w V_a T_w A_i%Environmental variables.
global  dx dy  %Grid variables
global   T_berg  %Physical constants


        %Deciding which velocity to use
        r_ind(1)=round((r(1)/dx))+1;
        r_ind(2)=round((r(2)/dy))+1;
        u_w_temp=[U_w(r_ind(1),r_ind(2)) V_w(r_ind(1),r_ind(2))];
        u_a_temp=[U_a(r_ind(1),r_ind(2)) V_a(r_ind(1),r_ind(2))];
        T_w_temp =T_w(r_ind(1),r_ind(2));
        A_i_temp =A_i(r_ind(1),r_ind(2));  %Sea ice area
 
       S_s=((2/3)*sqrt(norm(u_a_temp-u_w_temp))) +((1/10)*norm(u_a_temp-u_w_temp));%Sea State
       
       M_b=0.58*((norm(u_w_temp-u)).^(0.8)).*((T_w_temp-T_berg))./(norm(L).^(0.2))/(60*60*24); %Melt at the base    
       M_e=(1/12)*S_s*(1+cos(pi*(A_i_temp^(3))))*(T_w_temp+2)/(60*60*24) ; %Melt by wave erosion
       M_v=-((7.62*(10^(-3))*T_w_temp)+((1.29*(10^(-3))*(T_w_temp.^2))))/(60*60*24); %Melt by bouyant convection
        
