function [S,AR,qf,dS,qc]=geometry(i_panel,j_panel,c,b,Lambda,dx,dy,dxw,alpha,alphaV,panel_width)
%  
% Args:  
%   panel_width : dy
%   alphaV : twisting
%   alpha  : angle of attack
%   i_panel: no. of chordwise boxes
%   j_panel : no. of spanwise boxes

%%*******************************************************
% Original coding obtained from J. Katz, in fortran 77
% and found in :
% Low-Speed Aerodynamics,Joseph Katz et Allen Plotkin
% Second Edition,Cambridge Aerospace Series, 2001
% 
% translation in Matlab by Christophe Airiau, 2015
% free use of the matlab code, under J. Katz agreement.
%%*******************************************************
ib1=i_panel+1;
jb1=j_panel+1;
if length(alphaV) >= ib1
    Sin=sin(alphaV(1:ib1));
    Cos=cos(alphaV(1:ib1));
else
    fprintf('%s \n','Problem with dimension of alpha in geometry');
    return;
end

ctg1=tan(pi/2-Lambda(1));
ctg2=tan(pi/2-Lambda(2));
ctip=c+b*(ctg2-ctg1);
S=b*(c+ctip)/2;
AR=2*b*b/S;

%
%     wing fixed vortices location   ( qf(i,j,(x,y,z))...)
%
bj=0;
for j=1:jb1
    if j > 1 
        bj=bj+panel_width(j-1);
    end
    z1=0;
    dc1=bj*ctg1;
    dc2=bj*ctg2;
    dx1=(c+dc2-dc1)/i_panel;
    %     dc1=leading edge x,   dc2=trailing edge x
    for  i=1:i_panel
        qf(i,j,1)=dc1+dx1*(i-0.75);             % x coordinate
        qf(i,j,2)=bj;                           % y coordinate
        qf(i,j,3)=z1-0.25*dx1*Sin(i);           % z coordinate
        z1=z1-dx1*Sin(i);
    end

    %     the following lines are due to wake distance from trailing edge
    qf(ib1,j,1)=c+dc2+dxw;
    qf(ib1,j,2)=qf(i_panel,j,2);
    qf(ib1,j,3)=z1-dxw*Sin(i_panel);
end
%
%     wing colocation points
%
for j=1:j_panel
    z1=0;
    bj=qf(1,j,2)+panel_width(j)/2;
    dc1=bj*ctg1;
    dc2=bj*ctg2;
    dbx1=(c+dc2-dc1)/i_panel;
    for i=1:i_panel
        qc(i,j,1)=dc1+dx1*(i-0.25);
        qc(i,j,2)=bj;
        qc(i,j,3)=z1-0.75*dx1*Sin(i);
        z1=z1-dx1*Sin(i);
        dS(i,j)=dx1*panel_width(j);       % area of each panel
    end
end
%
%     rotation of wing points due to alpha
%
Sin1=-sin(alpha);
Cos1=cos(alpha);
for i=1:ib1
    for j=1:jb1
        qf1=qf(i,j,1);
        qf(i,j,1)=qf1*Cos1-qf(i,j,3)*Sin1;
        qf(i,j,3)=qf1*Sin1+qf(i,j,3)*Cos1;
        if (i ~= ib1) && (j < jb1)
            qc1=qc(i,j,1);
            qc(i,j,1)=qc1*Cos1-qc(i,j,3)*Sin1;
            qc(i,j,3)=qc1*Sin1+qc(i,j,3)*Cos1;
        end
    end
end

end
