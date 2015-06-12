function [W]= Kernal(norm_r_dist,L_eff,B1)
%d=B1*(L_eff/2);
d=(L_eff/2);
%d=B1.*sqrt(7./10).*(L_eff./2); 
%Note, d= sqrt(7/10)*(L_eff/2) is the min d needed to make A<1

R=norm_r_dist./d;

%W=(15./(7.*pi.*(d.^2))).*((((2/3)-(R.^2)+(0.5.*(R.^3))).*(R<1)) +(((1./6).*((2-R).^3)).*(R>1 & R< 2)));


%Gaussian
%W= (1./(pi.*(d.^2))).*exp(-((R.^2)));

W=(1./(pi.*(d.^2))).*(R<1);