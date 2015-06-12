function [P_ia P_ia_times_ui]= relative_motion_damping(berg_count,u,u_pred,r_n,u_n,L_eff,M)

global  spring_forces_on
global Num_bergs p_ia q_ia c_ed Ks_A  dt_max


tol=0.1;

P_ia=[0 0;0 0];
I=[1 0;0 1]; %Identity matrix
P_ia_times_ui=[0 0];


p_ia_eff=p_ia;



if spring_forces_on==1
    
    for count=1:Num_bergs   %Effect of each other berg on berg_count
        if count~=berg_count
            factor=0.8;  %Allows icebergs to damp eachother before they feel the spring force. 
            %r_dist=-[(r_n(count,1)-r_n(berg_count,1)) (r_n(count,2)-r_n(berg_count,2))];
            r_dist=-factor*[(r_n(count,1)-r_n(berg_count,1)) (r_n(count,2)-r_n(berg_count,2))];
            norm_r_dist=norm(r_dist);%sqrt(  ((r_n(:,1)-r_n(berg_count,1)).^2) +  ((r_n(:,2)-r_n(berg_count,2)).^2));
            
            L_crit=(L_eff(count,1)+L_eff(berg_count,1))./2;
            norm_r_dist(find(norm_r_dist==0))=tol;
            r_dist_hat=r_dist./norm_r_dist;
            P=r_dist_hat'*r_dist_hat;
            
            
            if norm_r_dist<L_crit
                A_damp=min(overlap_area(L_eff(count,1)./2,L_eff(berg_count,1)./2,norm_r_dist),min((pi*((L_eff(count,1)/2).^2)),(2*pi*((L_eff(berg_count,1)./2).^2))));
                
                %This idea was to increase damping for trapped bergs - has been abandoned.
                %              [A_damp trapped_in_a_berg]=overlap_area(L_eff(count,1)./2,L_eff(berg_count,1)./2,norm_r_dist);
                %               trapped_in_a_berg=0;
                %               T_power=1;  %Coeffitient that dictakes how quickly you move from the regular damping to the fully overlapped damping.
                %               %T tells you whether you are in an overlap state, T=min(T_o,T_min)/T_min
                %               [A_o trapped_in_a_berg]=overlap_area(L_eff(count,1)./2,L_eff(berg_count,1)./2,norm_r_dist);
                %                A_min=min((pi*((L_eff(count,1)/2).^2)),(pi*((L_eff(berg_count,1)./2).^2)));
                %                T=(A_o./A_min)^T_power;
                %                if (A_o/A_min)<0.5
                %                end
                %               c_t= (1-T)+(c_ed*T);
                %               c_t= (1*(~trapped_in_a_berg))+(c_ed*(trapped_in_a_berg));
                %               c_t=2*sqrt(Ks_A); %Critical_value.
                
             %   c_t=0;
             turn_on_damp_limit=0;
             if turn_on_damp_limit==1
                if p_ia~=0
                     p_ia=min((0.8/(dt_max*(0.5*(norm(u_n(count,:)-u)+ norm(u_n(count,:)-u_pred))))),10^(10)) ;
                end

                if q_ia~=0
                     q_ia=min((0.8/(dt_max*(0.5*(norm(u_n(count,:)-u)+ norm(u_n(count,:)-u_pred))))),10^(10)) ;
                end
             end
                
                
                %  rho_thick=M(berg_count,1)./(pi*((L_eff(berg_count,1)/2).^2));
                rho_thick=min(M(berg_count,1)./(pi*((L_eff(berg_count,1)/2).^2)),M(count,1)./(pi*((L_eff(count,1)/2).^2)));
                
                %Damping on the radial component of velocity.
                p_ia_eff=p_ia.*A_damp.*rho_thick;
                P_ia=P_ia+(p_ia_eff.*(0.5*(norm(P*((u_n(count,:)-u)'))+norm(P*((u_n(count,:)-u_pred)')))).*P);
                P_ia_times_ui=P_ia_times_ui+(p_ia_eff.*(0.5*(norm(P*((u_n(count,:)-u)'))+norm(P*((u_n(count,:)-u_pred)')))).*(P*((u_n(count,:))'))');
                
                %Damping on the tangental component of velocity.
                q_ia_eff=q_ia.*A_damp.*rho_thick;
                P_ia=P_ia+(q_ia_eff.*(0.5*(norm((I-P)*((u_n(count,:)-u)'))+norm((I-P)*((u_n(count,:)-u_pred)')))).*(I-P));
                P_ia_times_ui=P_ia_times_ui+(q_ia_eff.*(0.5*(norm((I-P)*((u_n(count,:)-u)'))+norm((I-P)*((u_n(count,:)-u_pred)')))).*((I-P)*((u_n(count,:))'))');
                
                
            end
            
        end
    end
end

%C_ia=0;
%C_ia_times_ui=0;

end