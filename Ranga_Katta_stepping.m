
M=(rho_b.*Th.*l.*w.*N) ;
if time_count==1
    a_n=((f0)*[u_n(:,2) -u_n(:,1)]);
end

%u_star=u_n+((f0*dt/2)*[u_n(:,2) -u_n(:,1)]);
u_star=u_n+((dt/2)*a_n);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Runge Kutta - Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X1=r_n;
[F_field DL_eff]=Force_field(X1,u_n,L,L_eff,M);
for berg_count=1:Num_bergs
   % A1(berg_count,:)=accel(berg_count,X1,u_n,u_star(berg_count,:),dt/2,Th(berg_count,1),l(berg_count,1),w(berg_count,1),N(berg_count,1),M,near_a_wall(berg_count),L_eff,Ranga_not_Velet);
    [A1(berg_count,:) a_1(berg_count,:) ]=accel(berg_count,X1,u_n,u_star(berg_count,:),dt/2,Th(berg_count,1),l(berg_count,1),w(berg_count,1),N(berg_count,1),M,near_a_wall(berg_count),L_eff,Ranga_not_Velet);
end
%l1=dt*(u_n);
V1=(u_n);
DL1=DL_eff+(dt.*(C1.*L.*((alpha_max-(L_eff./L))./(alpha_max-1))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Runge Kutta - Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X2=r_n+(V1*dt/2);
[F_field DL_eff]=Force_field(X2,u_n+(dt*A1/2),L,(L_eff+DL1/2),M);
for berg_count=1:Num_bergs
    %A2(berg_count,:)=accel(berg_count,X2,u_n,u_star(berg_count,:),dt/2,Th(berg_count,1),l(berg_count,1),w(berg_count,1),N(berg_count,1),M,near_a_wall(berg_count),L_eff,Ranga_not_Velet);
    [A2(berg_count,:)  a_2(berg_count,:) ]=accel(berg_count,X2,u_n,u_star(berg_count,:),dt/2,Th(berg_count,1),l(berg_count,1),w(berg_count,1),N(berg_count,1),M,near_a_wall(berg_count),L_eff,Ranga_not_Velet);
end
%l2=dt*(u_star+(dt*(A1/2)));
V2=(u_star+(dt*(A1/2)));
DL2=DL_eff+(dt*(C1.*(L).*(((alpha_max-(L_eff+DL1)./L))/(alpha_max-1))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Runge Kutta - Part 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X3=r_n+(dt*V2/2);
[F_field DL_eff]=Force_field(X3,u_n+(dt*A2/2),L,L_eff+(DL2/2),M);
for berg_count=1:Num_bergs
    %A3(berg_count,:)= accel(berg_count,X3,u_n,u_star(berg_count,:),dt,Th(berg_count,1),l(berg_count,1),w(berg_count,1),N(berg_count,1),M,near_a_wall(berg_count),L_eff,Ranga_not_Velet);
    [A3(berg_count,:) a_3(berg_count,:)]= accel(berg_count,X3,u_n,u_star(berg_count,:),dt,Th(berg_count,1),l(berg_count,1),w(berg_count,1),N(berg_count,1),M,near_a_wall(berg_count),L_eff,Ranga_not_Velet);
end
%l3=dt*(u_star+(dt*A2/2));
V3=(u_star+(dt*A2/2));
DL3=DL_eff+dt*(C1*(L).*(((alpha_max-(L_eff+DL2)./L))/(alpha_max-1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Runge Kutta - Part 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X4=r_n+(dt*V3);
[F_field DL_eff]=Force_field(X4,u_n+(dt*A3),L,L_eff+DL3,M);
for berg_count=1:Num_bergs  
    %A4(berg_count,:)= accel(berg_count,X4,u_n,u_star(berg_count,:),dt,Th(berg_count,1),l(berg_count,1),w(berg_count,1),N(berg_count,1),M,near_a_wall(berg_count),L_eff,Ranga_not_Velet);
    [A4(berg_count,:) a_4(berg_count,:)]= accel(berg_count,X4,u_n,u_star(berg_count,:),dt,Th(berg_count,1),l(berg_count,1),w(berg_count,1),N(berg_count,1),M,near_a_wall(berg_count),L_eff,Ranga_not_Velet);
end
V4=(u_star+(dt*A3));
DL4=DL_eff+dt*(C1*(L).*(((alpha_max-(L_eff+DL3)./L))/(alpha_max-1)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Updating u and r %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_np1=u_star+(dt*((A1+(2*A2)+(2*A3)+A4))/6);   %Explicit part of Coriolis has been included in u_star.
r_np1=r_n+(dt*(V1+(2*V2)+(2*V3)+V4)/6);
%L_eff=L_eff+((DL1+(2*DL2)+(2*DL3)+DL4)/6);
L_eff=max(L_eff+((DL1+(2*DL2)+(2*DL3)+DL4)/6),L);  %Applies the constraint L_eff<= L.


%a_np1=((f0)*[u_np1(:,2) -u_np1(:,1)]);
a_np1=(a_1+(2*a_2)+(2*a_3)+a_4)/6;

%u_star=u_n+((dt/2)*a_n);

a_n=a_np1;
