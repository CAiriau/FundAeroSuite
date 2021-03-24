function [] =Influence_Coefficients(Point_i)
%     program no. 11: influence coeff. of a rectilinear source/doublet panel
%       ----------------------------------------------------------------------
%
%     this program calculates the influence of a rectilinear panel at an
%     arbitrary point. (program by lindsey browne, 1988).


%%*******************************************************
% Original coding obtained from J. Katz, in fortran 77
% and found in :
% Low-Speed Aerodynamics,Joseph Katz et Allen Plotkin
% Second Edition,Cambridge Aerospace Series, 2001
% 
% translation in Matlab by Christophe Airiau, 2015
% free use of the matlab code, under J. Katz agreement.
%%*******************************************************

% point of interest  : P_i=(x_i,y_i,z_i)

if nargin == 0
    Point_i=[2 3 1];
end

epsilon=1.e-06;
pnds=1.0;    % p N ds ?
%    input doublet and source strengths
doublet=1.0/(4.0*pi);
sigma=1.0/(4.0*pi);

%     square/flat panel
%     input coordinates
PointPanel=zeros(5,3);              % A_i points
PointPanel(1,:)=[-0.5 -0.5 0];
PointPanel(2,:)=[ 0.5 -0.5 0];
PointPanel(3,:)=[ 0.5  0.5 0];
PointPanel(4,:)=[-0.5  0.5 0];
PointPanel(5,:)=PointPanel(1,:);
%     mid-point at (0,0,0)
Mid_point=[ 0 0 0];
fprintf('Mid point         : %12.6f %12.6f %12.6f \n',Mid_point);

rj31 = 0.0;
cjk1 = 0.0;
V_source=zeros(1,3);                % velocity due to source
V_doublet=zeros(1,3);               % velovity due to doublet
V=zeros(1,3);
Point=Point_i-Mid_point;            % frame centered on the Mid point O
fprintf('point of interest : %12.6f %12.6f %12.6f \n',Point);
distance=sqrt(sum(Point.^2));       % = d(OP_i)
for i=1:3
    Pnl(i)=0.25*sum(PointPanel(1:4,i));  % center point of the panel C
end
fprintf('centered point    : %12.6f %12.6f %12.6f \n',Pnl);
Pn=Point-Pnl;                    % \vec CO 
norm_Pn=sqrt(sum(Pn.^2));

d1=PointPanel(3,:)-PointPanel(1,:);    % d1=A_1A_3
d2=PointPanel(4,:)-PointPanel(2,:);    % d2=A_2A_4
cr=cross(d1,d2);                       %A_1A_3 x A_2A_4 : unit normal to the panel plane
norm_cr=sqrt(sum(cr.^2));
area=0.5*norm_cr;
cn=cr/norm_cr;                         % n =  normalized A_1A_3 x A_2A_4
pnn=sum(cn.*Pn);
tcm=0.5*(PointPanel(3,:)+PointPanel(4,:)) -Pnl;  % mid(A_3A_4)-C
norm_tcm=sqrt(sum(tcm.^2));                      % distance  mid(A_3A_4)-C
cm=tcm/norm_tcm;                                 % m = unit direction of A_3A_4
cl=cross(cm,cn);                                 % m x n = l third uni direction

for j=1:4
    k=j+1;
    a(1:3)=Point(1:3)-PointPanel(j,1:3);        % A_iP_i
    b(1:3)=Point(1:3)-PointPanel(k,1:3);        % A_kP_i
    s=PointPanel(k,:)-PointPanel(j,:);
    norm_a=sqrt(sum(a.^2));
    norm_b=sqrt(sum(b.^2));
    norm_s=sqrt(sum(s.^2));
    %   source contribution  
    sm=dot(s,cm);
    sl=dot(s,cl);
    am=dot(a,cm);
    al=dot(a,cl);
    bm=dot(b,cm);
    All=am*sl-al*sm;
    d_tmp=norm_a+norm_b-norm_s;
    if (d_tmp > 0) && (norm_s > 0)
        rj3=log((norm_a+norm_b+norm_s)/d_tmp)/norm_s  ;
    else
        rj3=0;
    end
    pa=Pn(3)^2*sl + All*am;  
    pb=pa - All*sm;  
    rnum=sm*Pn(3)*(norm_b*pa - norm_a*pb);  
    dnom=pa*pb + Pn(3)^2*norm_a*norm_b*sm^2;
    if abs(Pn(3)) < epsilon
        de=0;
    else
        if rnum ~= 0
            de=atan2(rnum,dnom);
        else
            de=0;
        end
    end
    rj31 = rj31 - sigma*All*rj3 ; 
    cjk1 = cjk1 -doublet*de;
    V_source=V_source+sigma*(rj3*(sm*cl-sl*cm)+de*cn);
%   doublet contribution  
    vmod=(norm_a+norm_b)/(norm_a*norm_b*(norm_a*norm_b +dot(a,b) ));
    V_doublet=V_doublet + doublet*vmod*cross(a,b);
end

%   limiting cases
dtt=2*pi;
if  distance> 0.0 
    pnds=Pn(3)^2/distance;
end
if (pnds < epsilon) && (distance > eps)
    dtt=Pn(3)*area/sqrt(distance)/distance;
end
if abs(dtt) < abs(cjk1)
    cjk1=dtt;
end
if  distance < epsilon^2 
    cjk1=-2*pi;
end

%   total
cjk = cjk1;
bjk = rj31 - Pn(3)*cjk1;
V=V_doublet+V_source;

norm_V_source=sqrt(sum(V_source.^2));
norm_V_doublet=sqrt(sum(V_doublet.^2));
norm_V=sqrt(sum(V.^2));

fprintf('area of panel      = %12.8f \n',area);
fprintf('source (potential) = %12.8f \n',bjk);
fprintf('source (velocity)  : \n');
fprintf('Vx_source          = %12.8f \n',V_source(1));
fprintf('Vy_source          = %12.8f \n',V_source(2));
fprintf('Vz_source          = %12.8f \n',V_source(3));
fprintf('V_source           = %12.8f \n',norm_V_source);
fprintf('doublet (potential): %12.8f \n', cjk);
fprintf('doublet (velocity) : \n');
fprintf('Vx_doublet         = %12.8f \n',V_doublet(1));
fprintf('Vy_doublet         = %12.8f \n',V_doublet(2));
fprintf('Vz_doublet         = %12.8f \n',V_doublet(3));
fprintf('V_doublet          = %12.8f \n',norm_V_doublet);
fprintf('total velocity     : \n');
fprintf('Vx                 = %12.8f \n',V(1));
fprintf('Vy                 = %12.8f \n',V(2));
fprintf('Vz                 = %12.8f \n',V(3));
fprintf('V                  = %12.8f \n',norm_V);

end


