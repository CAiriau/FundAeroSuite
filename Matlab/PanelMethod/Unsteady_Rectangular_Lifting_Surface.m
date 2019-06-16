function []=Unsteady_Rectangular_Lifting_Surface()
%     program no. 14: unsteady rectangular lifting surface (vlm)
%       
%
%     this is a 3-d linear code for rectangular planforms (with ground effect) 
%     in unsteady motion using the vortex lattice method (by joe katz, 1975).

%                       modes of operation
%
%     1. steady state        : set dt=dx/Vinf*i_panel*10,  and nsteps=5
%     2. sudden acceleration : set dt=dx/Vinf/4. and nsteps= up to 50
%     3. haeving oscillations: bh=heaving ampl. om= freq. 
%     4. pitch oscillations  : omega=frequency, teta=momemtary angle
%     5. for computational economy the parameter iw might be used
%     note; induced drag calculation increases computation time and
%     can be disconnected.
%%*******************************************************
% Original coding obtained from J. Katz, in fortran 77
% and found in :
% Low-Speed Aerodynamics,Joseph Katz et Allen Plotkin
% Second Edition,Cambridge Aerospace Series, 2001
% 
% translation in Matlab by Christophe Airiau, 2015
% free use of the matlab code, under J. Katz agreement.
%%*******************************************************
%
%
%     input data
%

global i_panel j_panel
%
% ==========
% input data
% ==========
%
line='====================================';
i_panel=4;           % number of panel in x direction
j_panel=13;          % number of panel in x direction
nsteps=5;           % number of time steps
nw=5;                % the number of (timewise) deforming wake elements.
rho=1.;              % density
bh=0;
omega=0;
Vinf=50;               % velocity m/s
c=1;                 % root chord
b=6;                 % semi-span
dx=c/i_panel;        % panel size in x
dy=b/j_panel;        % panel size in y
ch=10*c;             % ground clearance
alpha_1=5;           %
alpha_0=0;           %

alpha=pi/180*(alpha_1+alpha_0);
alf(1:i_panel)=0;
alf(i_panel+1)=0;     %  alf(i_panel+1) is required only for qf(i,j,k) calculation in geometry
% je ne sais pas ce que c'est alf ...

% alamda(i) are swep back angles. (alamda < 90, sweep backward)
Lambda(1)=90.*pi/180;
Lambda(2)=Lambda(1);
%dt=dx/Vinf/4.0;
dt=dx/Vinf*i_panel*10;
t=-dt;                  % time in seconds
dxw=0.3*Vinf*dt;
panel_width=dy*ones(1,j_panel);
%
%     initializations
%
ww=zeros(1,i_panel*j_panel);
dlt=zeros(i_panel,j_panel);
vortic0=zeros(i_panel,j_panel);
vort1=zeros(i_panel,j_panel);
Gamma=ones(i_panel,j_panel);    %  gama(i,j)=1 is required for influence matrix calculations
%
%     calculation of colocation points 
%
% dS : surface of a panel (i,j)
[S,AR,qf,dS,qc]=geometry(i_panel,j_panel,c,b,Lambda,dx,dy,dxw,0,zeros(1,i_panel+1),panel_width);
% geometry calculates wing colocation points qc,and vortex tips qf
%  maybe alpha instead of 0 .... to investigate
fprintf('Angle of attack            = %12.8f °\n',alpha*180/pi);
fprintf('wing semi-span             = %12.8f  \n',b);
fprintf('wing mean chord            = %12.8f  \n',c);
fprintf('wing aera                  = %12.8f  \n',S);
fprintf('wing aspect ratio          = %12.8f  \n',AR);
fprintf('free stream velocity       = %12.8f  \n',Vinf);
fprintf('height above ground        = %12.8f  \n',ch);
fprintf('number of panel in x       = %4i     \n',i_panel);
fprintf('number of panel in span(y) = %4i     \n',j_panel);
fprintf('size of dS  %d %d  \n', size(dS));
alf=alf*pi/180;
ib1=i_panel+1; jb1=j_panel+1;
%     program start

for it=1:nsteps     % time loop
    fprintf('%s\n',line);
    fprintf('time index                 = %4i     \n',it);
    fprintf('%s\n',line);
    t=t+dt; time(it)=t;
    sx=-Vinf*t;
    dsx=-Vinf;
    if ch >= 100
        ch1=0;
    else
        ch1=ch;
    end
    sz=bh*sin(omega*t)+ch1;
    dsz=bh*omega*cos(omega*t);
    % dsx=dsx/dt   dsz=dsz/dt
    teta=0.0;
    omega=0.0;             % ????? why here, depends on option ?
    Vinf=-cos(teta)*dsx-sin(teta)*dsz;
    sn1=sin(teta); cs1=cos(teta);
    wt=sn1*dsx-cs1*dsz;
    for i=1:i_panel
      sno(i)=sin(alpha+alf(i));
      cso(i)=cos(alpha+alf(i));
    end
    %
    %     ===========================
    %     vortex wake shedding points
    %     ===========================
    %
    for  j=1:jb1
        qw(it,j,1)=qf(ib1,j,1)*cs1-qf(ib1,j,3)*sn1+sx;
        qw(it,j,2)=qf(ib1,j,2);
        qw(it,j,3)=qf(ib1,j,1)*sn1+qf(ib1,j,3)*cs1+sz;
    end
    %
    %     ========================
    %     aerodynamic calculations
    %     ========================
    %
    %     influence coefficients calculation
    k=0;
    fprintf('Downwash ');
    for i=1:i_panel
        for j=1:j_panel
            fprintf('.');
            Sign=0.0;
            k=k+1;
            if it == 1
            %     matrix coefficients calculation occours only once for the
            %     time-fixed-geometry wing.
                w11=0;
                [a1,V]=wing0([qc(i,j,1) qc(i,j,2) qc(i,j,3)],Gamma,qf,sno,cso,Sign);
                l=0;
                for  i1=1:i_panel
                    for  j1=1:j_panel
                        l=l+1;
                        %     a(k,l) - is the normal velocity component due to a unit vortex lattice.
                        a(k,l)=a1(i1,j1);
                    end
                end
                % add influence of wing's other half 
                [a1,V]=wing0([qc(i,j,1) -qc(i,j,2) qc(i,j,3)],Gamma,qf,sno,cso,Sign);
                l=0;
                for i1=1:i_panel
                    for j1=1:j_panel
                        l=l+1;
                        a(k,l)=a(k,l)+a1(i1,j1);
                    end
                end
                if  ch <= 100.0 
                    % add influence of mirror image (due to ground)
                    Sign=10.0;
                    xx1=qc(i,j,1)*cs1-qc(i,j,3)*sn1+sx;
                    zz1=qc(i,j,1)*sn1+qc(i,j,3)*cs1+sz;
                    xx2=(xx1-sx)*cs1+(-zz1-sz)*sn1;
                    zz2=-(xx1-sx)*sn1+(-zz1-sz)*cs1;
                    [a1,V]=wing0([ xx2 qc(i,j,2) zz2],Gamma,qf,sno,cso,Sign);
                    l=0;
                    for i1=1:i_panel
                        for j1=1:j_panel
                            l=l+1;
                            a(k,l)=a(k,l)+a1(i1,j1);
                        end
                     end
                    %  add mirror image influence of wing's other half
                    [a1,V]=wing0([ xx2 -qc(i,j,2) zz2],Gamma,qf,sno,cso,Sign);
                    l=0;
                    for i1=1:i_panel
                        for j1=1:j_panel
                            l=l+1;
                            a(k,l)=a(k,l)+a1(i1,j1);
                        end
                    end
                    Sign=0;
                end
            else
                %     calculate wake influence
                xx1=qc(i,j,1)*cs1-qc(i,j,3)*sn1+sx;
                zz1=qc(i,j,1)*sn1+qc(i,j,3)*cs1+sz;
                [V]=wake([xx1 qc(i,j,2) zz1],it,qw,Gamma);
                [V1]=wake([xx1 -qc(i,j,2) zz1],it,qw,Gamma);
                if ch <= 100
                    [V2]=wake([xx1 qc(i,j,2) -zz1],it,qw,Gamma);
                    [V3]=wake([xx1 -qc(i,j,2) -zz1],it,qw,Gamma);
                else
                   V2=zeros(1,3);V3=zeros(1,3);
                end
                %    wake induced velocity is given in inertial frame
                V(1)=V(1)+V1(1)+V2(1)+V3(1);
                V(3)=V(3)+V1(3)-V2(3)-V3(3);
                u11= V(1)*cs1+V(3)*sn1;
                w11=-V(1)*sn1+V(3)*cs1;
                %     ww(k) is the prependicular component of wake influence to wing
                ww(k)=u11*sno(i)+w11*cso(i);
            end 

            %
            % calculate wing geometrical downwash
            %
            dw(k)=-Vinf*sno(i)+qc(i,j,1)*omega-wt;
            %   for general motion dw(k)=-Vinf*sin(alfa)+omega*x                    
            wts(i,j)=w11;
            %  w11 - is positive since the latest unsteady wake element is
            %  included in subroutine wing
        end % for j
    end % for i
    fprintf('/\n');

    %
    % solution of the problem:  dw(i)=ww(i)+a(i,j)*gama(i)              
    %
    k1=i_panel*j_panel;
    gama1(1:k1)=dw(1:k1)-ww(1:k1);
    Gamma1=a\gama1';
%     wing vortex lattice listing
    l=0;
    for i=1:i_panel 
        for j=1:j_panel
            l=l+1;
            Gamma(i,j)=Gamma1(l);
        end
    end
    %
    %     wake shedding
    %
    for j=1:j_panel
        %     latest wake elements listing
        vortic0(it,j)=Gamma(i_panel,j); vortic0(it+1,j)=0;
    end

    %
    %     ===========================
    %     wake distorsion calculation
    %     ===========================
    %
    iw=1;
    if it > 1 
        if it >= nw
            iw=it-nw+1;
            % nw is the number of (timewise) deforming wake elements.
        end
        i1=it-1;
        fprintf('Wake     ');
        for  i=iw:i1
            for j=1:jb1
                fprintf('.');
                [V]=veloce(qw(i,j,1:3),it,ch,Sign,alpha,sx,sz,qf,qw,sno,cso,vortic0);
                uvw(i,j,1:3)=V*dt;
            end
        end
        fprintf('/\n');
        %
        qw(iw:i1,1:jb1,1:3)=qw(iw:i1,1:jb1,1:3)+uvw(iw:i1,1:jb1,1:3);
    end
    %
    %     ==================
    %     Aerodynamic forces
    %     ==================
    %
    Lift=0;
    Drag=0;
    Pitching=0;
    Dynamic_Pressure=0.5*rho*Vinf^2;
    Gamma_total=0;
    fprintf('Forces   ');
    for j=1:j_panel
        fprintf('.');
        sigma=0; sigma1=0; dly(j)=0;
        for i=1:i_panel
            if i == 1
                gamaij=Gamma(i,j);              %  = gamma(x,y)
            else
                gamaij=Gamma(i,j)-Gamma(i-1,j);  % = d gamma / dx
            end
            dxm=(qf(i,j,1)+qf(i,j+1,1))/2;
            %     dxm is vortex distance from leading edge                       
            sigma1=(0.5*gamaij+sigma)*dx;
            sigma=Gamma(i,j);
            dfdt=(sigma1-dlt(i,j))/dt;
            %  dfdt  is the velocity potential time derivative
            dlt(i,j)=sigma1;
            dLift(i,j)=rho*(Vinf*gamaij+dfdt)*panel_width(j)*cso(i);
            % induced drag calculation 
            [V1]=wingl(qc(i,j,1:3),Gamma,qf);
            [V2]=wingl([ qc(i,j,1) -qc(i,j,2) qc(i,j,2)],Gamma,qf);
            if ch <= 100
                xx1=qc(i,j,1)*cs1-qc(i,j,3)*sn1+sx;
                zz1=qc(i,j,1)*sn1+qc(i,j,3)*cs1+sz;
                xx2=(xx1-sx)*cs1+(-zz1-sz)*sn1;
                zz2=-(xx1-sx)*sn1+(-zz1-sz)*cs1;
                [V3]=wingl([ xx2 qc(i,j,2) zz2],Gamma,qf);
                [V4]=wingl([ xx2 -qc(i,j,2) zz2],Gamma,qf);
            else
               V3(3)=0; V4(3)=0;
            end
            w8=V1(3)+V2(3)-V3(3)-V4(3);
            % add influence of mirror image (ground).
            cts=-(wts(i,j)+w8)/Vinf;
            dd1=rho*panel_width(j)*dfdt*sno(i);
            dd2=rho*panel_width(j)*Vinf*gamaij*cts;
            dDrag(i,j)=dd1+dd2;
            dp(i,j)=dLift(i,j)/dS(i,j)/Dynamic_Pressure;
            dly(j)=dly(j)+dLift(i,j);
            Lift=Lift+dLift(i,j);
            Drag=Drag+dDrag(i,j);
            Pitching=Pitching+dLift(i,j)*dxm;
            Gamma_total=Gamma_total+gamaij*panel_width(j);
        end
    end
    fprintf('/\n');
    Cl=Lift/(Dynamic_Pressure*S);
    Cd=Drag/(Dynamic_Pressure*S);
    Cm=Pitching/(Dynamic_Pressure*S*c);
    Cl_ref=2.*pi*alpha/(1.+2./AR);
    if abs(Cl_ref) < 1.e-20
        Cl_ref=Cl;
    end
    Cl_adim=Cl/Cl_ref;
    Gamma_ref=0.5*Vinf*S * Cl_ref;
    Gamma_adim=Gamma_total/Gamma_ref;
    %
    %     ======
    %     output
    %     ======
    %
    %        output
    %
    fprintf('Lift                       = %12.8f \n',Lift);
    fprintf('Drag                       = %12.8f \n',Drag);
    fprintf('Cl                         = %12.8f \n',Cl);
    fprintf('Cl_ref                     = %12.8f \n',Cl_ref);
    fprintf('Cd                         = %12.8f \n',Cd);
    fprintf('Cm                         = %12.8f \n',Cm);
    fprintf('f=Cl/Cd                    = %12.8f \n',Cl/Cd);
    fprintf('Cl_adim                    = %12.8f \n',Cl_adim);
    fprintf('Gamma total                = %12.8f \n',Gamma_total);
    fprintf('Gamma ref                  = %12.8f \n',Gamma_ref);
    fprintf('Gamma adim                 = %12.8f \n',Gamma_adim);
    fprintf('t                          = %12.8f \n',t);
    fprintf('sx                         = %12.8f \n',sx);
    fprintf('sz                         = %12.8f \n',sz);
    fprintf('Vinf                       = %12.8f \n',Vinf);
    fprintf('Theta                      = %12.8f \n',teta);
    fprintf('Omega                      = %12.8f \n',omega);

    for j=1:j_panel
        Delta_gammaj(2:i_panel,j)=diff(Gamma(:,j));
    end
    dly=dly./panel_width;
    if mod(it,5) == 0
        file=['output_time_' int2str(it) '.out'];
        fid=fopen(file,'w');
        for j=1:j_panel
            fprintf(fid,'j  = %d, dLy = %12.5f \n',j,dly(j));
            fprintf(fid,'delta p = %12.5f,%12.5f,%12.5f, %12.5f \n',dp(1:i_panel,j));
            fprintf(fid,'Gamma   = %12.5f,%12.5f,%12.5f, %12.5f \n',Gamma(1,j),Delta_gammaj(2:i_panel,j));
        end
        fprintf('\n Wake elements \n ');
        figure();
        N=size(vortic0);I=1:1:N(1);
        plot(I,vortic0);
        title('Gamma of each element');xlabel('time index');
        Mqw=size(qw);Mv=size(vortic0);
        for i=1:it
            fprintf(fid,'time                       = %12.8f \n',i*dt);
            fprintf(fid,'Gamma of each element\n');
            for k=1:Mv(2)
                fprintf(fid,'Gamma(%3i,%3i)   = %12.5f \n',i,k,vortic0(i,k));
            end
            for k=1:Mqw(2)
                fprintf(fid,'qw(%3i,%3i,:)= %12.8f %12.8f %12.8f \n',i,k,qw(i,k,:));
            end
        end
        fclose(fid);
    end
end   % temporal loop

end  % function

