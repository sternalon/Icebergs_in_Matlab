function [Force DL_total]= Force_field(r_n,u_n,L,L_eff,M)
global W_max Num_bergs Bond dx Land x y
global Ks C_damp Ks_land C_damp_land Ks_bond C_damp_bond C_p Ks_A rho_b dt_max
global alpha_max gamma mu B1 Nx Ny spring_forces_on use_land_force Ks_trapped PE_n


F_field_x=zeros(Num_bergs,1);
F_field_y=zeros(Num_bergs,1);
DL_total=zeros(Num_bergs,1);  %Total change to the size of all icebergs
PE_n=zeros(Num_bergs,1);

%P=(((L_eff-L).*(L_eff>L))./((alpha_max.*L)-L_eff)).^gamma;  %Percent being used for plastic deformation
P=(((L_eff-L).*(L_eff>L))./((alpha_max.*L)-L)).^gamma;  %Percent being used for plastic deformation


F_const=100000;
%tol2=0.000000000001;
%Note - the 2 should not be here, just temporarily
%F_e_max=C_p./((((L_eff./L)-1).*(L_eff>L))+tol2);

%alpha=3; %alpha>=1
Conc_n=zeros(Num_bergs,1);
for count=1:Num_bergs
    R=sqrt(  ((r_n(:,1)-r_n(count,1)).^2) +  ((r_n(:,2)-r_n(count,2)).^2));
    Conc_n(count)=sum(pi.*(L(:,1)./2).^(2).*Kernal(R,L_eff,B1));
end


tol=0.1;
%alpha=3;%0.01;%0.5;%%0.8;  %Gives you the value of C for

% %L_crit=W_max;
%
% %Calculating the concentration
% C=zeros(Num_bergs,1);
% for count=1:Num_bergs
%   % r_dist=[(r_n(:,1)-r_n(count,1)) (r_n(:,2)-r_n(count,2))];
%    norm_r_dist=sqrt(  ((r_n(:,1)-r_n(count,1)).^2) +  ((r_n(:,2)-r_n(count,2)).^2));
%    C(count,1)=sum(pi.*(L(:,1)./2).^(2).*Kernal(norm_r_dist,L,alpha)); %.*(norm_r_dist<L)
% end
% eps=0.00001;
% C(find(C>1))=1-eps;
% Gamma=1000*1;
% C0=0.8;
% %P=Gamma.*((C-C0)./(1-C));  %Calculating the pressure, from Gutfriaind
% C_star=50;
% P_star=1000000;%10^(6);
% P=P_star.*exp(-C_star.*(1-C));  %From Hibler
%
% P(find(P<0))=0;
% max(C)
%
% for count=1:Num_bergs
%    r_dist=[(r_n(:,1)-r_n(count,1)) (r_n(:,2)-r_n(count,2))];
%    norm_r_dist=sqrt(  ((r_n(:,1)-r_n(count,1)).^2) +  ((r_n(:,2)-r_n(count,2)).^2));
%    F_p_x(count)=sum(pi.*L(:,1).*L(:,1).*((P./(C.^2))+(P(count,1)./(C(count,1).^2))).*GradKernal(norm_r_dist,L,alpha).*(r_dist(:,1)))';
%    F_p_y(count)=sum(pi.*L(:,1).*L(:,1).*((P./(C.^2))+(P(count,1)./(C(count,1).^2))).*GradKernal(norm_r_dist,L,alpha).*(r_dist(:,2)))';
%    %F_p_x(count)=-sum(((P./(C.^2))+(P(count,1)./(C(count,1).^2))).*(2.*(r_dist(:,1))))';
%    %F_p_y(count)=-sum(((P./(C.^2))+(P(count,1)./(C(count,1).^2))).*(2.*(r_dist(:,2))))';
% end

if spring_forces_on==1
    
    %Only storing the forces at the relevant points:
    for count=1:Num_bergs
        r_dist=[(r_n(:,1)-r_n(count,1)) (r_n(:,2)-r_n(count,2))];
        norm_r_dist=sqrt(  ((r_n(:,1)-r_n(count,1)).^2) +  ((r_n(:,2)-r_n(count,2)).^2));
        
        %If particle length represents seperation distance
        %L_crit=(L(count,1)+L)./2;
        L_crit=(L_eff(count,1)+L_eff)./2;
        %If seperation distance is defined seperately
        %L_crit=W_max/3;
        
        no_self_force=ones(Num_bergs,1);
        no_self_force(count,1)=0;  %Make sure bergs do not influence them selves.
        norm_r_dist(find(norm_r_dist==0))=tol;
        
        %Elastic force for non-bonded particles  -Depending on distance
        %F_e=Ks*(L_crit-norm_r_dist).*(((norm_r_dist<L_crit)&~(Bond(count,:)'))).*(norm_r_dist>10^(-1));
        
        %Elastic force using overlapping areadnon-bonded particles
        
        %Creating a temporary spring constant.
        %       Ks_A_temp=Ks_A.*ones(Num_bergs,1);
        [A_o trapped_in_a_berg]=overlap_area(L_eff(count,1)./2,L_eff./2,norm_r_dist);
        A_o(count)=0;
        A_min=min((pi*((L_eff(count,1)/2).^2)),(pi*((L_eff./2).^2)));
        T=(A_o./A_min);
        N1=length(T);
        %     trapped_in_a_berg=T>.1;
        %     trapped_in_a_berg=max(T>.4).*ones(N1,1);
        %     Ks_A_temp=(Ks_A_temp.*(~trapped_in_a_berg))+(Ks_trapped.*trapped_in_a_berg);
        a=0.4; b=0.45;
        Ks_A_temp_coef= ((Ks_trapped.*(max(T)>b)) + (Ks_A.*(max(T)<a))+ (((max(T)<=b)&(max(T)>=a))*(((Ks_A-Ks_trapped)/(a-b)*(T-a))+Ks_A)));
        
        current_berg_trapped=(A_o==(pi*((L_eff(count,1)/2).^2)));
        
        %Case 1: we are considering forces on the small berg
        if max(current_berg_trapped)==1
            Ks_A_temp=Ks_A_temp_coef.*ones(N1,1);
        end
        %Case 2: we are considering forces on the bigger berg
        if max(current_berg_trapped)==0
            Ks_A_temp=(Ks_A.*(~current_berg_trapped))+(Ks_A_temp_coef.*current_berg_trapped);
        end
        %         Ks_A_temp=Ks_A.*ones(N1,1);
        %        Ks_A_temp=(Ks_A.*(~current_berg_trapped))+(Ks_A_temp_coef.*current_berg_trapped);
        
        
        rho_thick=M./(pi*((L_eff/2).^2));
        rho_thick_berg=M(count,1)./(pi*((L_eff(count,1)/2).^2));
        rho_thick_min=min(rho_thick_berg,rho_thick);
        %     F_e=rho_thick_min.*(Ks_A* overlap_area(L_eff(count,1)./2,L_eff./2,norm_r_dist).*(((norm_r_dist<L_crit)&~(Bond(count,:)'))).*(norm_r_dist>10^(-1)));
        F_e=rho_thick_min.*(Ks_A_temp.* A_o.*(((norm_r_dist<L_crit)&~(Bond(count,:)'))).*(norm_r_dist>10^(-1)));
        
        %Force using distance from center instead of over lap area - for analysis.
        K1=1;
        %F_e=rho_thick_min.*(K1.* (L_crit-norm_r_dist).*(((norm_r_dist<L_crit)&~(Bond(count,:)'))).*(norm_r_dist>10^(-1))).*no_self_force;
        PE_n=PE_n+0.5*(rho_thick_min.*(K1.* ((L_crit-norm_r_dist).^2).*(((norm_r_dist<L_crit)&~(Bond(count,:)'))).*(norm_r_dist>10^(-1)))).*no_self_force;
        
        
        
        %Elastic force for Bonded particles
        %  F_e=F_e+Ks_bond*(L_crit-norm_r_dist).*(((Bond(count,:)')));
        F_e=F_e+Ks_bond*(L_crit-norm_r_dist).*(((Bond(count,:)'))).*no_self_force;
        PE_n=PE_n+0.5*(0.5*Ks_bond*((L_crit-norm_r_dist).^2).*(((Bond(count,:)'))).*no_self_force);
        %note that we halve PE it because we are double counting each bond.
        
        
        
        
        
        
        F_s=F_e;
        
        %Turn on the dependence on concentration.
        %- Not working for very small bergs!!!
        %Only consider interaction force when Conc(r_n) is large
        %Defining interaction factor
        %q=(Conc_n).^(mu);   %First attempt
        %q=(max((Conc_n(count)+Conc_n-1),0).^(mu));    %Second attempt  -symetric. No interaction is the combined conc is less than 1
        %q=(((Conc_n(count)+Conc_n)/2).^(mu));   %Third attempt  -symetric.
        %F_s=F_s.*q;
        
        
        %Check if you are bigger than plastic limit. - Using > < instead of using an if statment
        %F_p=(abs(F_s)-F_e_max).*(F_s>F_e_max);
        %F_s=((-F_e_max).*(F_s>F_e_max))+(F_s.*(~(F_s>F_e_max)));
        
        F_p=F_s.*P;
        
        %still nedds to move the correct amount away from the wall.
        DL=-abs(F_p/Ks_A);  %Is it correct that it decreases the area of the circle, rather than something esle
        DL_total=DL_total+DL;
        
        %Adding Forces
        F_field_x = F_field_x +(F_s.*(r_dist(:,1)./norm_r_dist));
        F_field_y = F_field_y +(F_s.*(r_dist(:,2)./norm_r_dist));
        
    end
end
%end
DL(find(sum(Bond)>0))=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%[max(F_field_x) max(F_p_x)]

Force=zeros(Num_bergs,2);
Force(:,1)=F_field_x;%.*(M./(pi*((L_eff./2).^2)));
Force(:,2)=F_field_y;%.*(M./(pi*((L_eff./2).^2)));



