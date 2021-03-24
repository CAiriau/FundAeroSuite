function [] =Influence_Coefficients_old(x_i,y_i,z_i)
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
% point of interest  : x_i,y_i,z_i

epsilon=1.e-06;
pnds=1.0;    % p N ds ?
%    input doublet and source strengths
doublet=1.0/(4.0*pi);
sigma=1.0/(4.0*pi);

%     square/flat panel
%     input coordinates
x(1)=-.5;   y(1)=-.5;
x(2)=.5;    y(2)=-.5;
x(3)=.5;    y(3)=.5;
x(4)=-.5;   y(4)=.5;
z(1:4)=0;

x(5)=x(1);   y(5)=y(1);
z(5)=z(1);

%     mid-point at (0,0,0)
x0=0.0; y0=0.0; z0=0.0;

rj31 = 0.0;
cjk1 = 0.0;
Vx_source = 0.0; Vy_source = 0.0; Vz_source = 0.0;
Vx_doublet = 0.0;Vy_doublet = 0.0; Vz_doublet = 0.0;
Vx = 0.0;Vy = 0.0;Vz = 0.0;

px=x_i-x0; py=y_i-y0; pz=z_i-z0;
distance=sqrt(px*px+py*py+pz*pz);

pnlx=.25*sum(x(1:4)); pnly=.25*sum(y(1:4)); pnlz=.25*sum(z(1:4)); 
pnx=px-pnlx; pny=py-pnly; pnz=pz-pnlz;
pns=sqrt(pnx^2+pny^2+pnz^2);

d1x=x(3)-x(1); d1y=y(3)-y(1); d1z=z(3)-z(1);
d2x=x(4)-x(2); d2y=y(4)-y(2); d2z=z(4)-z(2);
crx=d1y*d2z-d2y*d1z; cry=d2x*d1z-d1x*d2z; crz=d1x*d2y-d2x*d1y;
crsq=sqrt(crx^2+cry^2+crz^2);
area=0.5*crsq;
cnx=crx/crsq; cny=cry/crsq; cnz=crz/crsq;
pnn=cnx*pnx+cny*pny+cnz*pnz;

tcmx=(x(3)+x(4))/2. - pnlx;
tcmy=(y(3)+y(4))/2. - pnly;
tcmz=(z(3)+z(4))/2. - pnlz;
tms=sqrt(tcmx^2+tcmy^2+tcmz^2);
cmx=((x(3)+x(4))/2. - pnlx)/tms;
cmy=((y(3)+y(4))/2. - pnly)/tms;
cmz=((z(3)+z(4))/2. - pnlz)/tms;
clx=cmy*cnz-cny*cmz;
cly=cnx*cmz-cmx*cnz;
clz=cmx*cny-cnx*cmy;

for j=1:4
    k=j+1;
    ax=px-x(j);     ay=py-y(j);     az=pz-z(j);
    bx=px-x(k);     by=py-y(k);     bz=pz-z(k);
    sx=x(k)-x(j);   sy=y(k)-y(j);   sz=z(k)-z(j);
    a=sqrt(ax*ax + ay*ay + az*az);
    b=sqrt(bx*bx + by*by + bz*bz);
    s=sqrt(sx*sx + sy*sy + sz*sz);
    %   source contribution  
    sm=sx*cmx+sy*cmy+sz*cmz;
    sl=sx*clx+sy*cly+sz*clz;
    am=ax*cmx+ay*cmy+az*cmz;
    al=ax*clx+ay*cly+az*clz;
    bm=bx*cmx+by*cmy+bz*cmz;
    All=am*sl-al*sm;
    if (a+b-s > 0) && (s > 0)
        rj3=log((a+b+s)/(a+b-s))/s  ;
    else
        rj3=0;
    end
    pa=pnz*pnz*sl + All*am;  
    pb=pa - All*sm;  
    rnum=sm*pnz*(b*pa - a*pb);  
    dnom=pa*pb + pnz*pnz*a*b*sm*sm;
    if abs(pnz) < epsilon
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
    Vx_source=Vx_source+sigma*(rj3*(sm*clx-sl*cmx)+de*cnx);
    Vy_source=Vy_source+sigma*(rj3*(sm*cly-sl*cmy)+de*cny);
    Vz_source=Vz_source+sigma*(rj3*(sm*clz-sl*cmz)+de*cnz);
%   doublet contribution  
    avbx=ay*bz - az*by;
    avby=az*bx - ax*bz;
    avbz=ax*by - ay*bx;
    adb=ax*bx + ay*by + az*bz;
    vmod=(a+b)/(a*b*(a*b + adb));
    Vx_doublet=Vx_doublet +doublet*vmod*avbx;
    Vy_doublet=Vy_doublet +doublet*vmod*avby;
    Vz_doublet=Vz_doublet +doublet*vmod*avbz;
end

%   limiting cases
dtt=2*pi;
if  distance> 0.0 
    pnds=pnz^2/distance;
end
if (pnds < epsilon) && (distance > eps)
    dtt=pnz*area/sqrt(distance)/distance;
end
if abs(dtt) < abs(cjk1)
    cjk1=dtt;
end
if  distance < epsilon^2 
    cjk1=-2*pi;
end

%   total
cjk = cjk1;
bjk = rj31 - pnz*cjk1;
Vx=Vx_doublet+Vx_source;
Vy=Vy_doublet+Vy_source;
Vz=Vz_doublet+Vz_source;

V_source=sqrt(Vx_source^2+Vy_source^2+Vz_source^2);
V_doublet=sqrt(Vx_doublet^2+Vy_doublet^2+Vz_doublet^2);
V=sqrt(Vx^2+Vy^2+Vz^2);

fprintf('area of panel      = %12.8f \n',area);
fprintf('source (potential) = %12.8f \n',bjk);
fprintf('source (velocity)  : \n');
fprintf('Vx_source          = %12.8f \n',Vx_source);
fprintf('Vy_source          = %12.8f \n',Vy_source);
fprintf('Vz_source          = %12.8f \n',Vz_source);
fprintf('V_source           = %12.8f \n',V_source);
fprintf('doublet (potential): %12.8f \n', cjk);
fprintf('doublet (velocity) : \n');
fprintf('Vx_doublet         = %12.8f \n',Vx_doublet);
fprintf('Vy_doublet         = %12.8f \n',Vy_doublet);
fprintf('Vz_doublet         = %12.8f \n',Vz_doublet);
fprintf('V_doublet          = %12.8f \n',V_doublet);
fprintf('total velocity     : \n');
fprintf('Vx                 = %12.8f \n',Vx);
fprintf('Vy                 = %12.8f \n',Vy);
fprintf('Vz                 = %12.8f \n',Vz);
fprintf('V                  = %12.8f \n',V);
 
end


