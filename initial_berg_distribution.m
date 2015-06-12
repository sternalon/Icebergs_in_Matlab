%clear all

%On off switches
consider_iceberg_density=0;  %Checks if new berg makes concentrations larger than 1
consider_tabular_berg_positions=0; %Checks if new tabular berg is too close to other bergs.
consider_iceberg_positions=0; %Checks if new berg is too close to other bergs

rng('default') %Initializes the random number generator
Num_groups=0;
Max_block_size=4;
block_size_x=ceil(Max_block_size*rand(Num_groups,1))+1;
%block_size_y=5;%2.*block_size_x;
block_size_y=ceil(Max_block_size*rand(Num_groups,1))+1;
Num_small_bergs=1;
%Defining iceberg initial position manually
just_one_line=1;
use_gridded_icebergs=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
Manual_bergs=0;
if use_gridded_icebergs==1
    MF=1.3/3; %Multiplication factor
    y_min= MF*W_max;%2*MF*W_max;%MF*W_max;
    y_max= Ly-y_min;%2*MF*W_max;%
    if just_one_line==1
        y_min=((y_min+y_max)./2); y_max=y_min;
    end
    x_min= 1.5*Lx/12;5*Lx/12;
    x_max=Lx-5*(MF*W_max);Lx-(MF*W_max);Lx/3+MF*W_max;
    y_grid=[y_min:MF*W_max:y_max];
    x_grid=[x_min:MF*W_max:x_max];
    Manual_bergs=(length(x_grid).*length(y_grid));
end

Num_bergs=sum(block_size_x.*block_size_y)+Num_small_bergs+Manual_bergs;
conc_threshold=1.005;
min_ind=1;  %Number of grid points between initial pos and land.

%Initializing all the iceberg vectors
r_n=zeros(Num_bergs,2);
N=ones(Num_bergs,1);%Ie: Number of bergs represented by a lagrangian point.Default is equal to 1.
L=zeros(Num_bergs,1); %Ie: the longer width
W=zeros(Num_bergs,1);%Ie: shorter width
Aspect_ratio=ones(Num_bergs,1); %Aspect ratio of little icebergs
l=zeros(Num_bergs,1); % Longer width of little icebergs inside a lagrangian point
w=zeros(Num_bergs,1); % Longer width of little icebergs inside a lagrangian point
L_eff=zeros(Num_bergs,1); %Ie: the longer width
Th=zeros(Num_bergs,1);
D=zeros(Num_bergs,1);
M=zeros(Num_bergs,1);
F=zeros(Num_bergs,1);
dr=zeros(Num_bergs,2);
lat=zeros(Num_bergs,1);
r_np1=zeros(Num_bergs,2);
u_np1=zeros(Num_bergs,2);
u_n=zeros(Num_bergs,2);
iceberg_death_flag=zeros(Num_bergs,1);
du=zeros(Num_bergs,2);
F_w=zeros(Num_bergs,2);
F_a=zeros(Num_bergs,2);
F_c=zeros(Num_bergs,2);
F_ss=zeros(Num_bergs,2);
F_r=zeros(Num_bergs,2);
F_i=zeros(Num_bergs,2);
a_1=zeros(Num_bergs,2);
a_2=zeros(Num_bergs,2);
a_3=zeros(Num_bergs,2);
k1=zeros(Num_bergs,2);
k2=zeros(Num_bergs,2);
k3=zeros(Num_bergs,2);
k4=zeros(Num_bergs,2);
l1=zeros(Num_bergs,2);
l2=zeros(Num_bergs,2);
l3=zeros(Num_bergs,2);
l4=zeros(Num_bergs,2);
DL1=zeros(Num_bergs,1);
DL2=zeros(Num_bergs,1);
DL3=zeros(Num_bergs,1);
DL4=zeros(Num_bergs,1);
last_good_position=zeros(Num_bergs,2);
m_e=zeros(Num_bergs,1);
m_v=zeros(Num_bergs,1);
m_b=zeros(Num_bergs,1);
Conc=zeros(Nx,Ny);
Bond=zeros(Num_bergs,Num_bergs);
L_Bond=zeros(Num_bergs,Num_bergs);
pos_group_x=zeros(Num_groups,1);
pos_group_y=zeros(Num_groups,1);
berg_number=zeros(Num_bergs,1); %Used to avoid bonding from seperate tabular bergs
a_n=zeros(Num_bergs,2);  %for Verlet Method
a_np1=zeros(Num_bergs,2); %for Verlet Method 
b_n=zeros(Num_bergs,2);  %for Verlet Method
b_np1=zeros(Num_bergs,2); %for Verlet Method 
u_star=zeros(Num_bergs,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Defining tabular icebergs %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if load_from_file==0

%Making icebergs begin as a block
berg_count=0;
if Num_groups>0
    berg_count=0;
    for k=1:Num_groups
        flag=0;
        berg_count_temp=berg_count;
        while flag==0
            flag=1;
            %pos_group_x(k)=Lx*rand(1,1);
            pos_group_x(k)=Lx/3+(3*Lx/8)*rand(1,1);
            pos_group_y(k)=(Ly/6)+((4*Ly/6)*rand(1,1));%Ly/6+(0*Lx/6)*rand(1,1);%Ly*rand(1,1);
            r_ind_1=round((pos_group_x(k)/dx))+1;
            r_ind_2=round((pos_group_y(k)/dy))+1;
            if Land(r_ind_1,r_ind_2)==1
                flag=0;
            else
                for k_past=1:k-1
                    if norm([pos_group_x(k) pos_group_y(k)]-[pos_group_x(k_past) pos_group_y(k_past)])< (Max_block_size*W_max);
                        if consider_tabular_berg_density==1
                            flag=0;
                        end
                    end
                end
            end
            
            
            for j=1:block_size_y(k)
                for i=1:block_size_x(k)
                    berg_count=berg_count+1;
                    %r_n(berg_count,1)=pos_group_x(k) +((i-1)*(2*W_max));
                    %r_n(berg_count,2)=pos_group_y(k) +((j-1)*(2*W_max));
                    r_n(berg_count,1)=pos_group_x(k) +((i-1)*(W_max));
                    r_n(berg_count,2)=pos_group_y(k) +((j-1)*(W_max));
                    berg_number(berg_count,1)=k;
                    % Setting iceberg size to be constant for now.
                    L(berg_count,1)=W_max; %Initial iceberg width
                    L_eff(berg_count,1)=W_max; %Initial iceberg width
                    W(berg_count,1)=W_max; %Initial iceberg width
                    Th(berg_count,1)=Th_max; %Initial iceberg thickness
                    [r_n(berg_count,1), r_n(berg_count,2)]=Enforce_boundary_conditions(r_n(berg_count,:));
                    
                    r_ind_1=round(r_n(berg_count,1)/dx)+1;
                    r_ind_2=round(r_n(berg_count,2)/dy)+1;
                    if Land(r_ind_1,r_ind_2)==1
                        flag=0;
                    end
                end
            end
            if flag==0
                berg_count=berg_count_temp;
            end
        end
    end
end

L_eff(find(N==1))=L(find(N==1));
%Calculating iceberg concentrations.
for i=1:Nx
    for j=1:Ny
        R=sqrt(  ((r_n(1:berg_count,1)-x(i)).^2) +  ((r_n(1:berg_count,2)-y(j)).^2));
        Conc(i,j)=sum(pi.*(L(1:berg_count,1)./2).^(2).*Kernal(R,L_eff(1:berg_count,1),B1));
        if Conc(i,j)>conc_threshold
            ['Conc > 1 at x= ' num2str(x(i)/1000) ', y= ' num2str(y(i)/1000)]
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Defining Singular icebergs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Defining position of single cell bergs
Bergs_so_far=berg_count;
k=Num_groups+1;
if Num_bergs-Manual_bergs>0
    %for berg_count=sum(block_size_x.*block_size_y)+1:Num_bergs
    for berg_count=Bergs_so_far+1:Num_bergs-Manual_bergs
        flag=0;
        
        % Define the initial size distribution of icebergs
        N(berg_count,1)=1;%randi(Max_Num-1,1)+1;
        %Macroscopic dimensions.
        L(berg_count,1)=W_max*rand(1,1);%(W_max/3)+(2*W_max/3.*rand(1,1)); %Initial iceberg width
        W(berg_count,1)=L(berg_count,1);%(W_max/2)+(W_max/2.*rand(1,1)); %Initial iceberg width
        Aspect_ratio(berg_count,1)=1;%0.3+(0.7*rand(1,1));
        Th(berg_count,1)=Th_max.*rand(1,1); %Initial iceberg thickness
        %L_eff(berg_count,1)=((alpha_max).*rand(1,1)).*L(berg_count,1);  %Initializing L_eff
        L_eff(berg_count,1)=1.0*L(berg_count,1);  %Initializing L_eff
        
        
        while flag==0
            flag=1;
            %r_n(berg_count,1)=Lx*rand(1,1);
            r_n(berg_count,1)=Lx/2+Lx/3*rand(1,1);%Lx/3+2*Lx/3*rand(1,1);;
            r_n(berg_count,2)=Ly/10+Lx/50*rand(1,1);;%+Ly/6*rand(1,1);%(Ly/6)+((4*Ly/6)*rand(1,1));%Ly/6*rand(1,1);
            
            %Checking if it is to near to land.
            r_up_ind=ceil((r_n(berg_count,:)/dx))+1;
            r_down_ind=floor((r_n(berg_count,:)/dx))+1;
            
            r_ind=round(r_n(berg_count,:)/dx)+1;
            for q= (r_ind(1)-min_ind):(r_ind(1)+min_ind)
                for j =(r_ind(2)-min_ind):(r_ind(2)+min_ind)
                    if q<1 | q>Nx | j<1  | j>Ny | Land(q,j)==1
                        flag=0;
                    end
                end
            end
            
            %         for q=[r_up_ind(1) r_up_ind(1) r_down_ind(1) r_down_ind(1)]
            %             for j=[r_up_ind(2) r_down_ind(2) r_up_ind(2) r_down_ind(2)]
            %                 if q<1 | j<1 |q>Nx | j>Ny
            %                     flag=0;
            %                 elseif Land(q,j)==1
            %                     flag=0;
            %                 end
            %             end
            %         end
            %         r_ind_1=round(r_n(berg_count,1)/dx)+1;
            %         r_ind_2=round(r_n(berg_count,2)/dy)+1;
            %         if Land(r_ind_1,r_ind_2)==1
            %             flag=0;
            %         else
            if consider_iceberg_positions==1
                for berg_past=1:berg_count-1
                    if (norm([r_n(berg_count,1) r_n(berg_count,2)]-[r_n(berg_past,1) r_n(berg_past,2)])< W_max);
                        flag=0;
                    end
                end
            end
            
            if consider_iceberg_density==1
                %Calculating iceberg concentrations.
                R=sqrt(  ((r_n(berg_count,1)-X).^2) +  ((r_n(berg_count,2)-Y).^2));
                Conc_temp=Conc+(pi.*(L(berg_count,1)./2).^(2).*Kernal(R,L_eff(berg_count,1),B1));
                if max(max(Conc_temp))>conc_threshold
                    flag=0;
                end
            end
        end
        if consider_iceberg_density==1
            Conc=Conc_temp;
        end
        
        
    end
    if consider_iceberg_density==1
        if max(max(Conc_temp))>1
            ['Conc > 1 inside intial_berg script'];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Defining Gridded icebergs %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if use_gridded_icebergs==1
    for i=1:length(x_grid)
        for j=1:length(y_grid)
            berg_count=berg_count+1;
            r_n(berg_count,:)=[x_grid(i) y_grid(j)];
        %    r_n(berg_count,:)=[x_grid(i) y_grid(j)]+(MF*W_max/5)*(rand(1,1)-0.5);%random offsets
            
            N(berg_count,1)=1;%randi(Max_Num-1,1)+1;
            L(berg_count,1)=W_max.*rand(1,1);%   (W_max/3)+(2*W_max/3.*rand(1,1)); %Initial iceberg width
            W(berg_count,1)=L(berg_count,1);%(W_max/2)+(W_max/2.*rand(1,1)); %Initial iceberg width
            Aspect_ratio(berg_count,1)=0.3+(0.7*rand(1,1));
            Th(berg_count,1)=Th_max;%.*rand(1,1); %Initial iceberg thickness
            L_eff(berg_count,1)=1.*L(berg_count,1);%+(((alpha_max-1).*rand(1,1))).*L(berg_count,1);  %Initializing L_eff
            %L_eff(berg_count,1)=*L(berg_count,1);  %Initializing L_eff
        end
    end
end



%Microscopic dimensions.
%Note that the aspect ratio of boned icebergs is automatically set to 1.
l=sqrt((pi*((L/2).^2)./(N.*Aspect_ratio)));
w=l.*Aspect_ratio;
%(N.*l.*w)-(pi.*(((L./2).^2)))


%Enforcing periodic boundary conditions
for berg_count=1:Num_bergs
    [r_n(berg_count,1), r_n(berg_count,2)]=Enforce_boundary_conditions(r_n(berg_count,:));
end
%
% %Defining initial bonds
for i=1:Num_bergs
    for j=1:Num_bergs
        %Bond(i,j)=(abs(norm(r_n(i,:)-r_n(j,:))-(2*W_max))<tol)*(berg_number(i,1)==berg_number(j,1));
        %L_Bond(i,j)=(abs(norm(r_n(i,:)-r_n(j,:))-(sqrt(2)*2*W_max))<tol)*(berg_number(i,1)==berg_number(j,1));
        Bond(i,j)=(abs(norm(r_n(i,:)-r_n(j,:))-(W_max))<tol)*(berg_number(i,1)==berg_number(j,1));
        L_Bond(i,j)=(abs(norm(r_n(i,:)-r_n(j,:))-(sqrt(2)*W_max))<tol)*(berg_number(i,1)==berg_number(j,1));
    end
end


%Using Short Bonds Only, getting rid of long bonds, including hexagical
for i=1:Num_bergs
    for j=1:Num_bergs
        if L_Bond(i,j)==1 & (mod((r_n(i,1)>r_n(j,1))+(r_n(i,2)>r_n(j,2)),2)==1  )
            Bond(i,j)=1;
        end
    end
end
L_Bond(:,:)=0;


if Jacobi_initial_bond_correction==1
    for k=1:100
        r_np1=r_n;
        jacobi_iteration;
        r_n=r_np1;
    end
    Jacobi_initial_bond_correction=0;
end


clear r_ind


% figure
% hold on
% axis([0 Lx/1000 0 Ly/1000])
% for berg_count=1:Num_bergs
%     hline(:,berg_count)=plot(r_n(berg_count,1)/1000,r_n(berg_count,2)/1000,'*r');
%     hline2(:,berg_count)=plot((r_n(berg_count,1)/1000)+(((L(berg_count,1)/2)/1000).*cos(circ_ind)),(r_n(berg_count,2)/1000)+(((L(berg_count,1)/2)/1000).*sin(circ_ind)),'g');
%     hline8(:,berg_count)=plot((r_n(berg_count,1)/1000)+(((L_eff(berg_count,1)/2)/1000).*cos(circ_ind)),(r_n(berg_count,2)/1000)+(((L_eff(berg_count,1)/2)/1000).*sin(circ_ind)),'m');
% end
