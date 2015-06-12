
r_np1=r_n+ (dt.*u_n)+ (((dt)^2/2)*a_n)+ (((dt)^2/2)*b_n);
%u_star=u_n+(0.5*dt*a_n)+((f0*dt/2)*[u_n(:,2) -u_n(:,1)]);
%u_star=u_n+((dt/2)*(a_n+(f0*[u_n(:,2) -u_n(:,1)])));

u_star=u_n+((dt/2)*(a_n));
DL=0;  %Note, the time stepping of DL has not been included in the Verlet scheme.

[F_field DL_eff]=Force_field(r_np1,u_n,L,L_eff+DL,M);
for berg_count=1:Num_bergs
    [Acc(berg_count,:) a_np1(berg_count,:) b_np1(berg_count,:) ]= accel(berg_count,r_np1,u_n,u_star(berg_count,:),dt,Th(berg_count,1),l(berg_count,1) ,w(berg_count,1),N(berg_count,1),M,near_a_wall(berg_count),L_eff,Ranga_not_Velet);
end
%u_np1=u_star+((dt/2)*(a_np1))+(dt*(b_np1));
u_np1=u_star+dt*Acc;

%Updating velocities and accelerations
a_n=a_np1;
b_n=b_np1;
r_n=r_np1;


