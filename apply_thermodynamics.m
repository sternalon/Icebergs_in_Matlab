
        for berg_count=1:Num_bergs
            %[m_b(berg_count,:) m_e(berg_count,:) m_v(berg_count,:)]=meltrate1(r_n(berg_count,:),u_n(berg_count,:),l(berg_count,1));  %Note, this depends on W in one of the terms
            m_b(berg_count,:) =0;
            m_e(berg_count,:)=0;
            m_v(berg_count,:)=0;
        end
        
        %Only melting the perimeter of tabular icebergs:
        l=l-(((m_e+m_v).*dt).*((6-sum(Bond)')./6));
        w=w-(((m_e+m_v).*dt).*((6-sum(Bond)')./6));
        Th=Th-(m_b*dt);
        
        L=2*sqrt(((N.*l.*w)/pi));
        L_eff(find(N==1))=L(find(N==1));        
        
        % Calculate the dimension and mass of iceberg
        M=(rho_b*Th.*l.*w).*N;
        D=(rho_b/rho_w).*Th; %Draft of iceberg
        F=Th-D; %Freeboard of iceberg
        
        if Rolling_icebergs==1
            for berg_count=1:Num_bergs
                %Iceberg rolling, from Martin and Adcroft - taken from Bigg et al and Weeks
                %Only for non-tabular icebergs
                if sum(Bond(berg_count,:))==0  %Alon's idea to not allow tabular bergs to roll
                    if l(berg_count,1)<sqrt((0.92.*(D(berg_count,1).^2))+(58.32*D(berg_count,1)))
                        temp=w(berg_count,1);
                        w(berg_count,1)=Th(berg_count,1);
                        Th(berg_count,1)=temp;
                    end
                end
            end
        end