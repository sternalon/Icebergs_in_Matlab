function [r1 r2]= Enforce_boudary_condition(r)
   
global  Lx Ly  %Grid variables

if (r(1)>(2*Lx)) | (r(2)>(2*Ly))
   halt 
end

        if r(1)>Lx
            r(1)=r(1)-Lx;
        end
        if r(2)>Ly
            r(2)=r(2)-Ly;
        end
        
        if r(1)<0
            r(1)=r(1)+Lx;
        end
        if r(2)<0
            r(2)=r(2)+Ly;
        end
        r1=r(1);
        r2=r(2);
   
        if (r1<0) | (r1>Lx)
            1
            halt
        end
        if r2<0 |r2>Ly
            2
            halt
        end
