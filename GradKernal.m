function [DW_Dx]= GradKernal(norm_r_dist,L,alpha)


%Calculates DW/Dx   (not including the (xi-xk) or (yi-yk) term.

d=alpha.*sqrt(7./10).*(L./2);
R=norm_r_dist./d;
%size(d)
%size(exp(-((R.^2))))

%W= (1./(pi.*(d.^2))).*exp(-((R.^2)));


%GradW=(15/(7*pi*d.^2)).*((-2.*R)+((3/2).*(R.^2)).*(R<1)) +((-(1/2).*((2-R).^2)).*(R>1 & R<2));

DW_Dx= (15./(7*pi*(d.^3))).*(((-2)+((3/2).*(R)).*(R<1)) +((-(1./(2*(R+(R==0)))).*((2-R).^2)).*(R>1 & R<2)));

