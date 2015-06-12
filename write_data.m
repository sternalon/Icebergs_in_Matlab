if time_count==1
    write_count=1;
    %Initialize variables for saving
    Write_total=ceil(end_time/write_period)+2;
    berg_pos_x=zeros(Num_bergs,Write_total);
    berg_pos_y=zeros(Num_bergs,Write_total);
    berg_vel_x=zeros(Num_bergs,Write_total);
    berg_vel_y=zeros(Num_bergs,Write_total);
    berg_pos_L=zeros(Num_bergs,Write_total);
    berg_pos_L_eff=zeros(Num_bergs,Write_total);
    berg_pos_l=zeros(Num_bergs,Write_total);
    berg_pos_w=zeros(Num_bergs,Write_total);
    berg_pos_N=zeros(Num_bergs,Write_total);
    berg_pos_Th=zeros(Num_bergs,Write_total);
    berg_pos_time=zeros(Write_total,1);
    
    %Saving initial values
    berg_pos_x(:,write_count)=r_n(:,1);
    berg_pos_y(:,write_count)=r_n(:,2);
    berg_vel_x(:,write_count)=u_n(:,1);
    berg_vel_y(:,write_count)=u_n(:,2);
    berg_pos_L(:,write_count)=L(:,1);
    berg_pos_L_eff(:,write_count)=L_eff(:,1);
    berg_pos_l(:,write_count)=l(:,1);
    berg_pos_w(:,write_count)=w(:,1);
    berg_pos_N(:,write_count)=N(:,1);
    berg_pos_Th(:,write_count)=Th(:,1);
    berg_pos_time(write_count)=time;
end

if (mod(time,dt*round(write_period/dt)))==0
    ['Writing data at time ' num2str(time/3600/24) ' (days)']
    write_count=write_count+1;
    berg_pos_x(:,write_count)=r_n(:,1);
    berg_pos_y(:,write_count)=r_n(:,2);
    berg_vel_x(:,write_count)=u_n(:,1);
    berg_vel_y(:,write_count)=u_n(:,2);
    berg_pos_L(:,write_count)=L(:,1);
    berg_pos_L_eff(:,write_count)=L_eff(:,1);
    berg_pos_l(:,write_count)=l(:,1);
    berg_pos_w(:,write_count)=w(:,1);
    berg_pos_N(:,write_count)=N(:,1);
    berg_pos_Th(:,write_count)=Th(:,1);
    berg_pos_time(write_count)=time;
end

%Write data to file
%     if write_iceberg_data==1
%         %Creating file on the first iteration.
%         if time_count==1 %& berg_count==1
%             fid=fopen(['data/iceberg' num2str(berg_count) '.txt'],'w');
%             fprintf(fid,'%5s\n',['time (days)      x (km)          y (km)         L (m)             W (m)            Th (m)       u (m/s)          v (m/s)   ']);
%             fclose(fid);
%         end
%         if (mod(time,dt*round(write_period/dt)))==0
%             fid=fopen(['data/iceberg' num2str(berg_count) '.txt'],'a');
%             fprintf(fid,'%5f %15f %15f %15f %15f %15f %15f  %15f\n',[(time/(60*60*24)) r_n(berg_count,1)/1000 r_n(berg_count,2)/1000 L(berg_count,1) W(berg_count,1) Th(berg_count,1) u_n(berg_count,1) u_n(berg_count,2)]);
%             fclose(fid);
%         end
%     end