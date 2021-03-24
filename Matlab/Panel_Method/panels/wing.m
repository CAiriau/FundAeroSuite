function [a1,V] = wing(Pt,Gamma,qf,ds,Sign,OnOff,i_panel,j_panel,i1,j1)
%        calculates induced velocity at a point (x,y,z), due to vorticity
%        distribution gamma(i,j), of semi-configuration - in a wing fixed
%        coordinate system.


%%*******************************************************
% Original coding obtained from J. Katz, in fortran 77
% and found in :
% Low-Speed Aerodynamics,Joseph Katz et Allen Plotkin
% Second Edition,Cambridge Aerospace Series, 2001
% 
% translation in Matlab by Christophe Airiau, 2015
% free use of the matlab code, under J. Katz agreement.
%%*******************************************************
V=0;
ib1=i_panel+1;
for i=1:ib1
    for j=1:j_panel
    % i3 is wake vortex counter
    i3=i;
    if i==ib1 
        i3=i_panel;
    end
    vortic=Gamma(i3,j);
    [V2]=vortex(Pt,qf(i,j+1,1:3),qf(i+1,j+1,1:3),vortic);
    [V4]=vortex(Pt,qf(i+1,j,1:3),qf(i,j,1:3),vortic);
    if OnOff >= 0.1  
        [V1]=vortex(Pt,qf(i,j,1:3),qf(i,j+1,1:3),vortic);
        [V3]=vortex(Pt,qf(i+1,j+1,1:3),qf(i+1,j,1:3),vortic);
        V0=V2+V4+V1+V3;
    else
        V1=zeros(1,3); V3=zeros(1,3);
        V0=V2+V4;
    end
    %fprintf('size %s : %d %d \n','V1',size(V1) );
    %fprintf('%d %d %d ( %d %d ) \n',size(ds(i1,j1,1:3)),i1,j1);

    a1(i,j)=V0(1)*ds(i1,j1,1)+V0(2)*ds(i1,j1,2)+V0(3)*ds(i1,j1,3);
    if Sign >= 1.0
        a1(i,j)=V0(1)*ds(i1,j1,1)+V0(2)*ds(i1,j1,2)-V0(3)*ds(i1,j1,3);
    end
    if i==ib1
        a1(i_panel,j)=a1(i_panel,j)+a1(ib1,j);
    end
    V=V+V0;
    end
end

end

