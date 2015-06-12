function [A trapped_in_a_berg]= overlap_area(R1,R2,d)
%Gives overlap between two circles. If one circle is inside the other, it
%gives the area of the smaller cirle. If there is no overlap, it returns
%zero. The flag trapped in a berg is returned when one circle is trapped
%inside the other.
A=((R1.^(2)).*acos(((d.^(2))+(R1.^(2))-(R2.^(2)))./(2.*d.*R1)))+ ((R2.^(2)).*acos(((d.^(2))+(R2.^(2))-(R1.^(2)))./(2.*d.*R2)))- (0.5.*sqrt( (-d+R1+R2).*(d+R1-R2).*(d-R1+R2).*(d+R1+R2) ));
A=real(A).*(d<(R1+R2));
% 

trapped_in_a_berg=(d<=norm(R1-R2));
A=A.*(~trapped_in_a_berg)+(trapped_in_a_berg.*(min(pi.*(R1.^2),pi*(R2.^2))));  %Requirement that one iceberg is not inside the other.
% if (d<=norm(R1-R2))
%     A=min(pi*(R1^2),pi*(R2^2));  %Requirement that one iceberg is not inside the other.
%     trapped_in_a_berg=1;
% end

