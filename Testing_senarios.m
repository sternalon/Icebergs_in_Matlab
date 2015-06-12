%
% % %One particle testing
if one_particle_testing==1
r_n(1,1)=Lx/2;r_n(1,2)=Ly/3;
%r_n(2,1)=Lx/2+dx/1 ;r_n(2,2)=Ly/2+dx/10;
Th(1,1)=10;% Th(2,1)=.001;
L(1,1)=W_max;% L(2,1)=W_max;
W=L; L_eff=L;
l=sqrt((pi*((L/2).^2)./(N.*Aspect_ratio)));
w=l.*Aspect_ratio;
u_n(1,1)=-0.5;
end

% %Two particle testing
if two_particle_testing==1
r_n(1,1)=Lx/3;r_n(2,1)=Lx/3+0.5*W_max;2*dx;
r_n(1,2)=Ly/2 ;r_n(2,2)=Ly/2;%+dx/10;
Th(1,1)=10; Th(2,1)=10;
L(1,1)=W_max; L(2,1)=W_max;
W=L; L_eff=L;
l=sqrt((pi*((L/2).^2)./(N.*Aspect_ratio)));
w=l.*Aspect_ratio;
u_n(1,1)=0; u_n(2,1)=-0;
Turn_all_damping_off=1;
Bond(2,1)=1;Bond(1,2)=1;
end

% %Three particle testing
if three_particle_testing==1
r_n(1,1)=Lx/2;r_n(2,1)=Lx/2+1*dx;r_n(3,1)=Lx/2+(3*W_max);
r_n(1,2)=Ly/2; r_n(2,2)=Ly/2; r_n(3,2)=Ly/2;%+dx/10;
Th(1,1)=10; Th(2,1)=10; Th(3,1)=10;
L(1,1)=3*W_max; L(2,1)=W_max/5; L(3,1)=3*W_max;
W=L; L_eff=L;
l=sqrt((pi*((L/2).^2)./(N.*Aspect_ratio)));
w=l.*Aspect_ratio;
u_n(1,1)=0.0; u_n(2,1)=-0.0; u_n(3,1)=-0.0;
Turn_all_damping_off=1;
end

%Three particle testing against wall.
if three_particle_testing2==1
r_n(1,1)=Lx/8;r_n(2,1)=Lx/8+1*dx;r_n(3,1)=Lx/8+3*dx;
r_n(1,2)=Ly/2; r_n(2,2)=Ly/2; r_n(3,2)=Ly/2;%+dx/10;
Th(1,1)=10; Th(2,1)=10; Th(3,1)=10;
L(1,1)=W_max; L(2,1)=W_max/2; L(3,1)=W_max;
W=L; L_eff=L;
l=sqrt((pi*((L/2).^2)./(N.*Aspect_ratio)));
w=l.*Aspect_ratio;
u_n(1,1)=-1.0; u_n(2,1)=-1.0; u_n(3,1)=-1.0;
end


% %THree particle testing_bonded
if three_particle_testing_bonded==1
r_n(1,1)=Lx/3;r_n(2,1)=Lx/3+0.5*W_max;2*dx; r_n(3,1)=Lx/3-0.5*W_max;
r_n(1,2)=Ly/2 ;r_n(2,2)=Ly/2;r_n(3,2)=Ly/2;
Th(1,1)=.00010; Th(2,1)=10; Th(3,1)=10;
L(1,1)=W_max; L(2,1)=W_max; L(3,1)=W_max;
W=L; L_eff=L;
l=sqrt((pi*((L/2).^2)./(N.*Aspect_ratio)));
w=l.*Aspect_ratio;
u_n(1,1)=0; u_n(2,1)=-0; u_n(3,1)=-0;
Turn_all_damping_off=0;
%Bond(2,1)=1;Bond(1,2)=1; Bond(1,3)=1;Bond(3,1)=1;
end