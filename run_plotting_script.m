
if plot_just_velocities==1
    time_hist(time_count)=time/(60*60*24);
    % speed_hist(time_count)=norm(u_n(berg_view_num,:));
    % speed_hist(time_count)=(u_n(berg_view_num,1));
  %  speed_hist(time_count)=mean(((u_n(:,1).^2)+(u_n(:,2).^2)));    
    speed_hist(time_count)=u_n(berg_view_num,1);    

    figure(2);hold on
    if Ranga_not_Velet==1
        plot(time/(60*60*24),speed_hist(time_count),'b*');
        plot(time_hist,speed_hist,'b');
    end
    if Ranga_not_Velet==0
        plot(time/(60*60*24),speed_hist(time_count),'r*');
        plot(time_hist,speed_hist,'r');
    end
    xlabel('time (day)')
end
if plot_spring_energy==1
    time_hist(time_count)=time/(60*60*24);
    if plot_the_total_energy==0
        %Note: PE of one particle does not really make sense.
        PE_hist(time_count)=(PE_n(berg_view_num,1));
        KE_hist(time_count)=0.5*((M(berg_view_num,1).*((u_np1(berg_view_num,1).^2)+(u_np1(berg_view_num,2).^2))));
    end
    if plot_the_total_energy==1
        PE_hist(time_count)=sum(PE_n);
        KE_hist(time_count)=0.5*(sum((M(:,1).*((u_np1(:,1).^2)+(u_np1(:,2).^2)))));
    end
    
    
    E_hist=PE_hist+KE_hist;
    
    figure(2);hold on
    plot(time_hist,E_hist,'r','linewidth',4);
    plot(time_hist,PE_hist,'b:');
    plot(time_hist,KE_hist,'k:');

    %plot(time/(60*60*24),E_hist(time_count),'b*');
    xlabel('time (day)')
end


%     if plot_land_history==1
%         for land_count=1:Num_moving_land
%            land_hist(find(x<moving_land_r_n(land_count,1) & moving_land_u_n(land_count,1)~=0),)=1;
%         end
%        pcolor(x,y,land_hist);shading('interp')
%     end


if plot_the_area==1
    if mod(time,dt*round(plot_period/dt))==0
        %    Calculating area on a grid
        if time_count>1 %& berg_count==1
            delete(hline6);
        end
        for i=1:Nx
            for j=1:Ny
                R=sqrt(  ((r_n(:,1)-x(i)).^2) +  ((r_n(:,2)-y(j)).^2));
                Conc(i,j)=sum(pi.*(L(:,1)./2).^(2).*Kernal(R,L_eff,B1));
            end
        end
        if max(max(Conc))>1
            time_count;
            ['Conc > 1  in loop!'];
        end
    end
    hline6=pcolor(x/1000,y/1000,Conc');shading('interp');caxis([0.0 1]);colorbar
end

if plot_land_each_time==1
    if mod(time,dt*round(plot_period/dt))==0
        if time_count>1 %& berg_count==1
            delete(hline10);
            delete(hline11);
        end
        
        hold on
        for i=1:Nx
            for j=1:Ny
                if Land(i,j)==1
                    hline10(i,j)=plot(x(i)/1000,y(j)/1000,'k*');
                    hline11(i,j)=plot((x(i)/1000)+(((dx/2)/1000).*cos(circ_ind)),(y(j)/1000)+(((dx/2)/1000).*sin(circ_ind)),'k','linewidth',4);
                    
                end
            end
        end
    end
end

if plot_each_time==1
    figure(1)
    for berg_count=1:Num_bergs
        %Plotting figures
        if mod(time,dt*round(plot_period/dt))==0
            if plot_each_time==1
                %Displaying the position
                % if time_count>1 %& berg_count==1
                delete(hline(:,berg_count));
                delete(hline2(:,berg_count));
                delete(hline8(:,berg_count));
                % end
                hline(:,berg_count)=plot(r_n(berg_count,1)/1000,r_n(berg_count,2)/1000,'*r');
                hline2(:,berg_count)=plot((r_n(berg_count,1)/1000)+(((L(berg_count,1)/2)/1000).*cos(circ_ind)),(r_n(berg_count,2)/1000)+(((L(berg_count,1)/2)/1000).*sin(circ_ind)),'g');
                hline8(:,berg_count)=plot((r_n(berg_count,1)/1000)+(((L_eff(berg_count,1)/2)/1000).*cos(circ_ind)),(r_n(berg_count,2)/1000)+(((L_eff(berg_count,1)/2)/1000).*sin(circ_ind)),'m');
                
                
                if plot_the_bonds==1
                    %Plotting the bonds
                    for pair_berg=1:berg_count
                        if (Bond(pair_berg,berg_count)==1)
                            if time_count>1 %& berg_count==1
                                
                                delete(hline4(:,berg_count,pair_berg));
                            end
                            if  norm(r_n(berg_count,:)-r_n(pair_berg,:))<(Lx/2)
                                hline4(:,berg_count,pair_berg)=plot([r_n(berg_count,1)/1000 r_n(pair_berg,1)/1000],[r_n(berg_count,2)/1000 r_n(pair_berg,2)/1000],'m','linewidth',3);
                            end
                        end
                        %                             if (L_Bond(pair_berg,berg_count)==1)
                        %                                 if time_count>1 %& berg_count==1
                        %                                     delete(hline5(:,berg_count,pair_berg));
                        %                                 end
                        %                                 if  norm(r_n(berg_count,:)-r_n(pair_berg,:))<(Lx/2)
                        %                                     hline5(:,berg_count,pair_berg)=plot([r_n(berg_count,1)/1000 r_n(pair_berg,1)/1000],[r_n(berg_count,2)/1000 r_n(pair_berg,2)/1000],'r','linewidth',3);
                        %                                 end
                        %                             end
                    end
                end
            end
            
            
            if berg_count==1
                if time_count>1 %& berg_count==1
                    delete(hline3);
                end
                %  hline3=text((((x_max+x_min)/2)/2)-(((x_max+x_min)/2)/20),(((y_max+y_min)/2))+(((y_max+y_min)/2)/20),['time = ' num2str(round(10*time/60/60/24)/10) ' days'],'fontsize',fontsize);
                hline3=text((((x_max+x_min)/2)),y_max+((y_max-y_min)/30),['time = ' num2str(round(10*time/60/60/24)/10) ' days'],'fontsize',fontsize);
            end
        end
    end
    
    if use_moving_land==1 & plot_movgin_land==1
        for land_count=1:Num_moving_land
            if mod(time,dt*round(plot_period/dt))==0
                %Displaying the position
                % if time_count>1 %& berg_count==1
                %delete(hline70(:,land_count));
                %delete(hline71(:,land_count));
                % end
                hline70(:,land_count)=plot(moving_land_r_n(land_count,1)/1000,moving_land_r_n(land_count,2)/1000,'*k');
                hline71(:,land_count)=plot((moving_land_r_n(land_count,1)/1000)+(((moving_land_L(land_count,1)/2)/1000).*cos(circ_ind)),(moving_land_r_n(land_count,2)/1000)+(((moving_land_L(land_count,1)/2)/1000).*sin(circ_ind)),'k');
                
            end
        end
    end
end

%Deciding if the icebeg is completely gone.
if thermodynamics_on==1
    if (L(berg_count,1)<0) | (W(berg_count,1)<0) | (Th(berg_count,1)<0)
        iceberg_death_flag(berg_count)=1;
        
        %Deletes it from the figure.
        delete(hline(:,berg_count));
        delete(hline2(:,berg_count));
    end
end
drawnow