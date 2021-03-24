function [time,Cl,Cd,Cl_t,gam1]=lumped_vortex()
%    program no. 13: sudden acceleration of a flat plate (lumped vortex)
%       -------------------------------------------------------------------
%
%     transient aerodynamic of a flat plate represented by a single lumped
%     vortex element (prepared as a solution of homework problems by joe
%     katz, 1987).

%%*******************************************************
% Original coding obtained from J. Katz, in fortran 77
% and found in :
% Low-Speed Aerodynamics,Joseph Katz et Allen Plotkin
% Second Edition,Cambridge Aerospace Series, 2001
% 
% translation in Matlab by Christophe Airiau, 2015
% free use of the matlab code, under J. Katz agreement.
%%*******************************************************
close all;
title='sudden accelaration of a flat plate';
line='********************************************';
fprintf('%s\n%s\n%s\n',line,title,line);
fprintf('%s \n',' it   x          Cd        Cl       Clt      \Gamma ');
format=(' %d  %12.5e    %12.5e    %12.5e    %12.5e    %12.5e  \n');

%
%  DATA
% 
n_vortex=50;
vortic=zeros(n_vortex,3);
%       vortic(i,1:2) =(x,z) of vortex
%       vortic(i,3)   = Gamma_i
nstep=20;       % number of time step
rho=1.;         % density
ut=50;          % velocity m/s
c=1;            % chord
alpha=5;        % incidence in degress
alpha_rad=alpha/180*pi;
Sin=sin(alpha_rad);Cos=cos(alpha_rad);
dt=c/(4*ut);    % time step
t=-dt;
dxw=0.3*ut*dt;   % distance in the wake

%     program start

for it=1:nstep
    t=t+dt;
    time(it)=t;
    %     path of origin (sx,sz)
    sx=-ut*t; sz=0;
    %     shedding of wake points
    vortic(it,1)= (c+dxw)*Cos+sx; 
    vortic(it,2)=-(c+dxw)*Sin+sz;
    %
    %     calculate momentary vortex strength of wing and wake vortices
    %
    a=-1/(pi*c);
    b=1/(2.0*pi*(c/4+dxw));
    rhs2=0;
    wwake=0;
    if it > 1 
        it1=it-1;
        %     calculate wake influence
        xx1= 0.75*c*Cos+sx;
        zz1=-0.75*c*Sin+sz;
        [u,w]=downwash(xx1,zz1,1,it1,vortic);
        wwake=u*Sin+w*Cos;
        %     calculation of rhs
        rhs2=rhs2-sum(vortic(1:it1,3));
    end

    rhs1=-ut*Sin-wwake;
    %     solution (based on algebraic solution of two equations for gammat 
    %     and the latest wake vortex strength vortic(it,3).
    vortic(it,3)=1/(b/a-1)*(rhs1/a-rhs2);
    gammat=rhs2-vortic(it,3);
    %
    %     wake  rollup
    %
    if it >= 1
        for i=1:it
            xx1= 0.25*c*Cos+sx;
            zz1=-0.25*c*Sin+sz;
            [u,w]=vor2D(vortic(i,1)-xx1,vortic(i,2)-zz1,gammat);
            [u1,w1]=downwash(vortic(i,1),vortic(i,2),1,it,vortic);
            u=u+u1; w=w+w1;
            uw(i,1)=vortic(i,1)+u*dt;
            uw(i,2)=vortic(i,2)+w*dt;
        end
        vortic(1:it,1:2)=uw(1:it,1:2);
    end
    %
    %    aerodynamic loads
    %
    if it==1 
        gamat1=0;
    end
    q=0.5*rho*ut^2;
    dgamdt=(gammat-gamat1)/dt;
    gamat1=gammat;
    %     calculate wake induced downwash
    xx1= 0.75*c*Cos+sx;
    zz1=-0.75*c*Sin+sz;
    [u,w]=downwash(xx1,zz1,1,it,vortic);
    ww=u*Sin+w*Cos;
    Lift=rho* (ut*gammat+dgamdt*c);
    Drag=rho*(-ww*gammat+dgamdt*c*Sin);
    Cl(it)=Lift/(q*c);
    Cd(it)=Drag/(q*c);
    %
    %     output
    %
    Cl_t(it)=Cl(it)/(2.0*pi*Sin);
    gam1(it)=gammat/(pi*ut*c*Sin);
    sx1(it)=sx-ut*dt;
    fprintf(format,it,sx1(i),Cd(i),Cl(i),gam1(i),Cl_t(i));
end
figure(1);
plot(time,Cl,'b-o',time,Cl_t,'r-');
legend('Cl','Clt');xlabel('time');
%title('Cl and Clt versus time');

figure(2);
plot(time,Cd);
%title('Cd versus time');
ylabel('Cd');xlabel('time');

figure(3);
plot(time,gam1);
ylabel('\Gamma');xlabel('time');
%title('\Gamma versus time');

