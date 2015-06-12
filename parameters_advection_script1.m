


global  U_w U_a U_i V_w V_a V_i eta_x eta_y Th_i T_w A_i%Environmental variables
global  dx dy dt %Grid variables
global  rho_w c_w c_dw rho_a c_a c_da f0 beta g rho_i c_iv T_berg rho_b  %Physical constants
global C_p alpha_max gamma mu B1 Ks_A c_ia c_wall p_ia q_ia c_ed Ks_trapped dt_max

 %Plastic constant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Basic parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt_max=1*3600; %maximum time step allowed
Lx=500*1000;%100*1000; %Length of the domain
Ly=500*1000;%120*1000 %Width of the domain
dx=25*1000;  %ocean grid spacing in x
dy=dx; %ocean grid spacing in y
dt=.1*3600; %time step
start_time=0;  %start time
end_time=250*24*3600;%10000000;  %end time
tol=0.000000001;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Planetary parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Omega=7.292*10^(-5); %Anular velocity.
R_earth=6371*(10^(3)); %Radius of the earth..
g=-9.81;
lat_ref=(60/180)*pi; % Reference latitude for using a beta plane

%Defining Coriolis parameter and beta parameter
f0=2*Omega*sin(lat_ref);
beta=2*Omega*cos(lat_ref)/R_earth;
beta=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Iceberg parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Num_bergs=9;  %Number of icebergs being used
W_max=20*1000;%*dx;  %Maximum width of icebegs (m)  - must be smaller than dx
Th_max=10;  %Maximum depth of icebegs
rho_b=900; %Iceberg density (kg/m^3), Hunke and Comeau 2011
T_berg=-4;
Max_Num=100; %Maximum number of icebergs represented by a point.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Iceberg interaction parameters%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Interaction Force constants
%Ks=1000;%10000; %Interaction Spring Force coeffitient
%C_damp=500; %Interaction Damping coeffitient
%This is the forward Euler minimum spring constant for a crit damped spring
Ks_crit_damp=0.25*((-3+sqrt(4+(dt^2)))/((dt^2)-1));
Ks_A=1*10^(-4);  %1*10^(-3) is the best so far.  Tests show that Ks_A=0.01 is too large with c_ia=0.01
%Ks_trapped=Ks_A;
Ks_trapped=min(Ks_A,Ks_crit_damp)%10^(-5);2/(dt_max^2);%Ks_A/100;
if Ks_A==0
   'Warning Ks_A=0'
end
if (Ks_A*(dt^2)/2)>1
   'Warning, spring constant above stability condition'
   (Ks_A*(dt^2)/2)
end
if Ks_trapped~=2/(dt_max^2)
   'Warning, Ks_trapped above CFL' 
   Ks_trapped./(2/(dt_max^2))
end
if Ks_trapped>Ks_crit_damp
       'Warning, Ks_trapped above CFL for critically damped spring.' 
end

c_wall=10*10^(-4);; %Implicit drag when overlapping with another iceberg or near a wall.
c_ia=1*10^(-4);%2*sqrt(Ks_A); %Implicit drag when overlapping with another iceberg, all all velocity.
p_ia=1*10^(-4);1*10^(-4);5*10^(-4);%1*10^(2); %Implicit drag when overlapping with another iceberg, on radial component
q_ia=1*10^(-4); %Implicit drag when overlapping with another iceberg, on tengental component

%c_ed=1*10^(9);%2*sqrt(Ks_A);%; %Enhanced damping for completely overlapping icebergs.

%Note, critical damping should have p_ia=c_ia=2*sqrt(Ks_A);
critical_damping_on=1;
if Ks_A~=0  & critical_damping_on
   % c_ia=2*sqrt(Ks_A) %Critical_value.
    p_ia=2*sqrt(Ks_A); %Critical_value.
    q_ia=2*sqrt(Ks_A); %Critical_value.
end
%p_ia=100;
%q_ia=0;
c_ia=0;

Ks_bond=100000; %Bonded Spring Force coeffitient
C_damp_bond=5000; %Bonded Damping coeffitient

Ks_land=0*10^(-3); %Land Spring Force coeffitient
%C_damp_land=0.00*D1; %Land Damping coeffitient


%Parameters used for iceberg spreading
alpha_max=4;  %This gives the maximum size of a cloud, L_max=alpha_max*L.
%C1=0.0000001;%.01;  %This constant gives the constant spreading of iceberg clouds (should perhaps depend on N - number of bergs in the cloud.
C1=0;%1/(5.*60*60*24);%.01;  %This constant gives the constant spreading of iceberg clouds (should perhaps depend on N - number of bergs in the cloud.

%C_p=10;  %This constant is not used for now. Perhaps delete it.
gamma=1;%0.2;  %Exponent for working out what percent of force to use on platisity.
mu=5;  % Exponent for working out whether to make particles repel eachother.  mu=0 turns this off.
B1=1; %Multiplicative factor used for estimating the Kernal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Ocean parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_w=1027;  %Ocean density
c_w=0.85;      %Ocean drag on iceberg coefitient on iceberg walls, Hunke and Comeau 2011
c_dw=5*10^(-4);%Ocean drag coefitient on iceberg base, Hunke and Comeau 2011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Atmosphere parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_a=1.2754; %Google
c_a=0.4;%Atmosphere drag coefitient on iceberg base, Hunke and Comeau 20
c_da=2.5*10^(-4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Sea ice parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho_i=900 ; %Sea ice density (kg/m^3), Hunke and Comeau 2011
c_iv=0.85;      %sea ice drage coefitient set equal to ocean drag, which is from Hunke and Comeau 2011, rather than Martin and Adcoft

Number_nearest_neighbours=10;  %How many of the closest neighbours to use when using the Jacobi contact correstion.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Plotting and writing parameters %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontsize=14;
plot_period=1*3600;% How often to plot output (in seconds)
write_period=1*24*3600;
filename='Test_new_forces';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Plotting and writing parameters %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_filename='New_test_curved_coast';