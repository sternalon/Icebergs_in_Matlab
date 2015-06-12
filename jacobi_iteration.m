%% Jacobi land corrections
if Jacobi_correct_land==1
    for i=1:Num_bergs
        
        r_up_ind=ceil((r_np1(i,:)/dx))+1;
        r_down_ind=floor((r_np1(i,:)/dx))+1;
        for q=[r_up_ind(1) r_up_ind(1) r_down_ind(1) r_down_ind(1)]
            for j=[r_up_ind(2) r_down_ind(2) r_up_ind(2) r_down_ind(2)]
                if ((q)<(Nx+1)) && ((q)>0) && ((j)<(Ny+1)) && ((j)>0)
                    
                    %         L_crit=(L_eff(i,1)/2+dx/2);
                    %         Num_points=ceil(L_crit/dx);
                    %         r_ind=round(r_n(i,:)/dx)+1;
                    %         for q=r_ind(1)-Num_points:r_ind(1)+Num_points
                    %             if (q<(Nx+1)) & (q>0)
                    %                 for j=r_ind(2)-Num_points:r_ind(2)+Num_points
                    %                     if (j<(Ny+1)) & (j>0)
                    %
                    
                    if Land(q,j)==1
                        
                        r_diff=(r_np1(i,:)-[x(q) y(j)]);  %Not set up for land near the edges
                        
                        %%Find out if you are in case 1,2 or 3.
                        flag_x=0;
                        flag_y=0;
                        
                        dir=sign((r_n(i,1)-x(q)));
                        test_num=2;
                        
                        %   r_dist=[(r_n(count,1)-x(q)) 0];
                        %   dir=r_dist(1)/norm(r_dist(1));
                        if ((q+dir)<(Nx+1)) && ((q+dir)>0)
                            if Land((q+dir),j)~=1
                                flag_x=1;
                            end
                        end
                        
                        dir=sign((r_n(i,2)- y(j)));
                        
                        % r_dist=[0 (r_n(count,2)- y(j))];
                        % dir=r_dist(2)/norm(r_dist(2));
                        if ((j+dir)<(Ny+1)) && ((j+dir)>0)
                            if Land(q,(j+dir))~=1
                                flag_y=1;
                            end
                        end
                        
                        if flag_x==flag_y
                            r_diff=[r_np1(i,1)-x(q) r_np1(i,2)-y(j)];  %Not set up for land near the edges
                            
                            
                        elseif flag_x==0
                            %  r_diff=(r_np1(i,:)-[x(q) y(j)]);
                            r_diff=[0 r_np1(i,2)-y(j)];
                        else
                            % r_diff=(r_np1(i,:)-[x(q) y(j)]);
                            r_diff=[r_np1(i,1)-x(q) 0];
                        end
                        
                        %     error=norm(r_diff)-((dx/2)+(L(i,1)/2));
                        %  error=norm(r_diff)-((dx)+(L(i,1)/2))
                        %  error=norm(r_diff)-((L(i,1)/2))
                        error=-((dx/2)-(norm(r_diff)-((L_eff(i,1)/2))));
                        
                        if error<0
                            r_np1(i,:)= r_np1(i,:)-((error).*r_diff./norm(r_diff));
                            %u_np1
                            %lowering the velocity
                            
                            P=[1 0;0 1]-((r_diff/norm(r_diff))'*(r_diff/norm(r_diff)));
                            u_np1(i,:)=(P*u_np1(i,:)')';  %Killing velocity perpendicular to the wall.
                            %                                u_np1(i,:)=u_np1(i,:)+ (u_np1(i,:).*r_diff./norm(r_diff))
                            %  r_diff./norm(r_diff)
                        end
                    end
                end
            end
        end
    end
end


%% Jacobi bond corrections

if Jacobi_bond_correction==1  |  Jacobi_initial_bond_correction==1
    for i=1:Num_bergs
        
        %Restoring the periodic boundaries
        if use_periodic_boundaries==1
            [r_np1(i,1), r_np1(i,2)]=Enforce_boundary_conditions(r_np1(i,:));
            [r_n(i,1), r_n(i,2)]=Enforce_boundary_conditions(r_n(i,:));
            [r_np1(i,1), r_np1(i,2)]=Enforce_boundary_conditions(r_np1(i,:));
            
            
            %Finding initial iceberg index
            r_ind(i,:)=round((r_np1(i,:)/dx))+1;
            if r_ind(i,1)==Nx+1
                r_ind(i,1)=1;
            end
            if r_ind(i,2)==Ny+1
                r_ind(i,2)=1;
            end
            if r_ind(i,1)==-1
                r_ind(i,1)=Nx;
            end
            if r_ind(i,2)==-1
                r_ind(i,2)=Ny;
            end
            
        end
        
        
        %Restoring bonds and %contact force
        
        %Only using nearest_neighbour limit
        Bond_list=[find(Bond(i,:)==1) find(L_Bond(i,:)==1)];
        Bond_list=sort(Bond_list); %So that it does not do short then long bonds
        if ~isempty(Bond_list)
            for p=1:length(Bond_list)
                j=Bond_list(p);
                % for j=1:i   %Old version
                
                if i>j  %This makes sure that we do not double count.
                    
                    if Bond(i,j)==1
                        %                             %Setting lenghts to be equal
                        %                             L(i,1)=0.5*(L(i,1)+L(j,1));
                        %                             L(j,1)=L(i,1);
                        %                             W(i,1)=0.5*(W(i,1)+W(j,1));
                        %                             W(j,1)=W(i,1);
                        %                             Th(i,1)=0.5*(Th(i,1)+Th(j,1));
                        %                             Th(j,1)=Th(i,1);
                        
                        %Sorting out short bonds
                        r_diff=(r_np1(i,:)-r_np1(j,:));
                        if norm(r_diff(1))>(Lx/2)
                            if r_np1(i,1)>r_np1(j,1)
                                r_np1(j,1)=r_np1(j,1)+Lx;
                            elseif r_np1(j,1)>r_np1(i,1)
                                r_np1(i,1)=r_np1(i,1)+Lx;
                            end
                            r_diff=(r_np1(i,:)-r_np1(j,:));
                        end
                        if norm(r_diff(2))>(Ly/2)
                            if r_np1(i,2)>r_np1(j,2)
                                r_np1(j,2)=r_np1(j,2)+Ly;
                            elseif r_np1(j,2)>r_np1(i,2)
                                r_np1(i,2)=r_np1(i,2)+Ly;
                            end
                            r_diff=(r_np1(i,:)-r_np1(j,:));
                        end
                        error=norm(r_diff)-(0.5*(L(i,1)+L(j,1)));
                        %error=norm(r_diff)-(W_max);
                        r_np1(i,:)= r_np1(i,:)-(error/2).*r_diff./norm(r_diff);
                        r_np1(j,:)= r_np1(j,:)+(error/2).*r_diff./norm(r_diff);
                        if use_periodic_boundaries==1
                            [r_np1(i,1), r_np1(i,2)]=Enforce_boundary_conditions(r_np1(i,:));
                            [r_np1(j,1), r_np1(j,2)]=Enforce_boundary_conditions(r_np1(j,:));
                        end
                    end
                    
                    if L_Bond(i,j)==1
                        %Setting lenghts to be equal
                        %                             L(i,1)=0.5*(L(i,1)+L(j,1));
                        %                             L(j,1)=L(i,1);
                        %                             W(i,1)=0.5*(W(i,1)+W(j,1));
                        %                             W(j,1)=W(i,1);
                        %                             Th(i,1)=0.5*(Th(i,1)+Th(j,1));
                        %                             Th(j,1)=Th(i,1);
                        
                        %Sorting out short bonds
                        r_diff=(r_np1(i,:)-r_np1(j,:));
                        if norm(r_diff(1))>(Lx/2)
                            if r_np1(i,1)>r_np1(j,1)
                                r_np1(j,1)=r_np1(j,1)+Lx;
                            elseif r_np1(j,1)>r_np1(i,1)
                                r_np1(i,1)=r_np1(i,1)+Lx;
                            end
                            r_diff=(r_np1(i,:)-r_np1(j,:));
                        end
                        if norm(r_diff(2))>(Ly/2)
                            if r_np1(i,2)>r_np1(j,2)
                                r_np1(j,2)=r_np1(j,2)+Ly;
                            elseif r_np1(j,2)>r_np1(i,2)
                                r_np1(i,2)=r_np1(i,2)+Ly;
                            end
                            r_diff=(r_np1(i,:)-r_np1(j,:));
                        end
                        %error=norm(r_diff)-(sqrt(2)*W_max);
                        error=norm(r_diff)-(sqrt(((L(i,1))^2)+((L(j,1)^2))));
                        r_np1(i,:)= r_np1(i,:)-(error/2).*r_diff./norm(r_diff);
                        r_np1(j,:)= r_np1(j,:)+(error/2).*r_diff./norm(r_diff);
                        if use_periodic_boundaries==1
                            [r_np1(i,1), r_np1(i,2)]=Enforce_boundary_conditions(r_np1(i,:));
                            [r_np1(j,1), r_np1(j,2)]=Enforce_boundary_conditions(r_np1(j,:));
                        end
                    end
                end
            end
        end
    end
end

%% Jacobi contact correction (iceberg colisions)

if Jacobi_contact_correction==1
    for i=1:Num_bergs
        R=sqrt(  ((r_n(:,1)-r_n(i,1)).^2) +  ((r_n(:,2)-r_n(i,2)).^2));
        [temp dist_order]=sort(R);
        for outer_count=1:min(Number_nearest_neighbours,Num_bergs)
            j=dist_order(outer_count);
            %Making icebergs stay away from eachother
            if (Bond(i,j)==0) & (L_Bond(i,j)==0)
                if i~=j
                    r_diff=(r_np1(i,:)-r_np1(j,:));
                    if use_periodic_boundaries==1
                        if norm(r_diff(1))>(Lx/2)
                            if r_np1(i,1)>r_np1(j,1)
                                r_np1(j,1)=r_np1(j,1)+Lx;
                            elseif r_np1(j,1)>r_np1(i,1)
                                r_np1(i,1)=r_np1(i,1)+Lx;
                            end
                            r_diff=(r_np1(i,:)-r_np1(j,:));
                        end
                        if norm(r_diff(2))>(Ly/2)
                            if r_np1(i,2)>r_np1(j,2)
                                r_np1(j,2)=r_np1(j,2)+Ly;
                            elseif r_np1(j,2)>r_np1(i,2)
                                r_np1(i,2)=r_np1(i,2)+Ly;
                            end
                            r_diff=(r_np1(i,:)-r_np1(j,:));
                        end
                    end
                    error=norm(r_diff)-((L_eff(i,1)/2)+(L_eff(j,1)/2));% Distance between bergs must be greater than the sum of their biggest side.
                    if error<0
                        r_np1(i,:)= r_np1(i,:)-(error/2).*r_diff./norm(r_diff);
                        r_np1(j,:)= r_np1(j,:)+(error/2).*r_diff./norm(r_diff);
                        
                        %Killing the velocities afere adjustment -
                        %completely inelastic collisions (for now).
                        a_n(i,:)=0;
                        a_n(j,:)=0;
                        b_n(i,:)=0;
                        b_n(j,:)=0;
                        u_np1(i,:)=0;
                        u_np1(j,:)=0;
                        
                        %                     %We can think later about only killing the one
                        %                     %component of the velocity (for example)
                        %                          P=[1 0;0 1]-((r_diff/norm(r_diff))'*(r_diff/norm(r_diff)));
                        %                          u_np1(i,:)=(P*u_np1(i,:)')';  %Killing velocity perpendicular to the wall.
                        %                          u_np1(j,:)=(P*u_np1(j,:)')';
                        
                        
                        
                    end
                    if use_periodic_boundaries==1
                        [r_np1(i,1), r_np1(i,2)]=Enforce_boundary_conditions(r_np1(i,:));
                        [r_np1(j,1), r_np1(j,2)]=Enforce_boundary_conditions(r_np1(j,:));
                    end
                end
            end
        end
    end
end

%% Jacobi correction around moving land boundaries

if use_moving_land==1
    for i=1:Num_bergs
        R=sqrt(  ((moving_land_r_n(:,1)-r_n(i,1)).^2) +  ((moving_land_r_n(:,2)-r_n(i,2)).^2));
        [temp dist_order]=sort(R);
        for outer_count=1:min(Number_nearest_neighbours,Num_bergs)
            j=dist_order(outer_count);
            %Making icebergs stay away from eachother
            r_diff=(r_np1(i,:)-moving_land_r_n(j,:));
            
            error=norm(r_diff)-((L_eff(i,1)/2)+(moving_land_L(j,1)/2));% Distance between bergs must be greater than the sum of their biggest side.
            if error<0
                r_np1(i,:)= r_np1(i,:)-(error).*r_diff./norm(r_diff);
                %  r_np1(j,:)= r_np1(j,:)+(error/2).*r_diff./norm(r_diff);
            end
        end
    end
end

