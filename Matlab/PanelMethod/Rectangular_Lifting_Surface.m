function []=Rectangular_Lifting_Surface()
%   3d-vlm code for simple wing planforms with ground effect (by joe katz, 1974).

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
% ==========
% input data
% ==========
%
i_panel=4;           % number of panel in x direction
j_panel=13;          % number of panel in x direction
x=[ 0 0 4 4 ];
b=13;               % semi span
Vinf=1.0;             % free stream velocity
alpha=5.0;          % angle of attack in degrees
ch=1000. ;           % height above ground
%        x(1) to x(4) are x-coordinates of the wing's four cornerpoints.
%        b - wing span, Vinf - free stream speed, b - wing span, 
%        ch - height above ground
%        constants
dxw=100.0*b;        % length of the wake
Gamma=ones(i_panel,j_panel);
%        Gamma(i,j)=1  is required for influence matrix calculations.
rho=1.;  % density
%   alpha=alpha1*pay/180;
%   sn1=sin(alpha)
%   cs1=cos(alpha)
%   ib1=ib+1
%   ib2=ib+2
%   jb1=jb+1
%
%        =============
%        wing geometry
%        =============
%
[S,c,b,AR,qf,ds,qc]=grid(alpha,b,ch,dxw,x,i_panel,j_panel);
fprintf('Angle of attack            = %12.8f °\n',alpha);
fprintf('wing semi-span             = %12.8f  \n',b);
fprintf('wing mean chord            = %12.8f  \n',c);
fprintf('wing aera                  = %12.8f  \n',S);
fprintf('wing aspect ratio          = %12.8f  \n',AR);
fprintf('free stream velocity       = %12.8f  \n',Vinf);
fprintf('height above ground        = %12.8f  \n',ch);
fprintf('number of panel in x       = %4i     \n',i_panel);
fprintf('number of panel in span(y) = %4i     \n',j_panel);
fprintf('size of ds  %d %d %d \n', size(ds));
%
%        ========================
%        aerodynamic calculations
%        ========================
%
%        influence coefficients calculation   
%

k=0;
for i=1:i_panel
    for j=1:j_panel
        Sign=0;
        k=k+1;
        [a1,V]=wing(qc(i,j,1:3),Gamma,qf,ds,Sign,1.0,i_panel,j_panel,i,j);
        l=0;
        for i1=1:i_panel
            for j1=1:j_panel
                l=l+1;
                % a(k,l) - is the normal velocity component due to a unit vortex
                %     lattice.
                a(k,l)=a1(i1,j1);
            end
        end
        % add influence of wing's other half 
        Pt_tmp=qc(i,j,1:3); Pt_tmp(2)=-Pt_tmp(2);
        [a1,V]=wing(Pt_tmp,Gamma,qf,ds,Sign,1.0,i_panel,j_panel,i,j);
        l=0;
        for i1=1:i_panel
            for j1=1:j_panel
                l=l+1;
                a(k,l)=a(k,l)+a1(i1,j1);
            end
        end
        if  ch <= 100.0 
            fprintf('ground effect');
            % add influence of mirror image (due to ground)
            Sign=10.0;
            Pt_tmp=qc(i,j,1:3); Pt_tmp(3)=-Pt_tmp(3);
            [a1,V]=wing(Pt_tmp,Gamma,qf,ds,Sign,1.0,i_panel,j_panel,i,j);
            l=0;
            for i1=1:i_panel
                for j1=1:j_panel
                    l=l+1;
                    a(k,l)=a(k,l)+a1(i1,j1);
                end
            end
            %  add mirror image influence of wing's other half
            Pt_tmp=-qc(i,j,1:3); Pt_tmp(1)=-Pt_tmp(1);
            [a1,V]=wing(Pt_tmp,Gamma,qf,ds,Sign,1.0,i_panel,j_panel,i,j);
            l=0;
            for i1=1:i_panel
                for j1=1:j_panel
                    l=l+1;
                    a(k,l)=a(k,l)+a1(i1,j1);
                end
            end
            Sign=0;
        end
        %
        %        calculate wing geometrical downwash
        %
        Uinf=[ Vinf 0  0];
        %        this is the general formulation for right hand side.
        n=ds(i,j,1:3);
        dw(k)=-dot(Uinf,reshape(n,1,3));
    end
end

%
%        solution of the problem:  dw(i)=a(i,j)*gamma(i)
%
Gamma1=a\dw';
l=0;
for i=1:i_panel 
    for j=1:j_panel
        l=l+1;
        Gamma(i,j)=Gamma1(l);
    end
end
%  to test : reshape(Gamma1,j_panel,i_panel)'
%
%        ==================
%        forces calculation
%        ==================
%
Lift=0;
Drag=0;
Pitching=0;
Dynamic_Pressure=0.5*rho*Vinf^2;
%dly=zeros(1,j_panel);

for  j=1:j_panel
    dly(j)=0;
    for i=1:i_panel
        if i==1   
            gamaij=Gamma(i,j);
        else
            gamaij=Gamma(i,j)-Gamma(i-1,j);
        end
        dym=qf(i,j+1,2)-qf(i,j,2);
        dl(i,j)=rho*Vinf*gamaij*dym;
        %        induced drag calculation 
        [a1,V1]=wing(qc(i,j,1:3),Gamma,qf,ds,Sign,0.0,i_panel,j_panel,i,j);
        Pt_tmp=qc(i,j,1:3); Pt_tmp(2)=-Pt_tmp(2);
        [a1,V2]=wing(Pt_tmp,Gamma,qf,ds,Sign,0.0,i_panel,j_panel,i,j);
        if ch <= 100.0 
            Pt_tmp=qc(i,j,1:3); Pt_tmp(3)=-Pt_tmp(3);
            [a1,V3]=wing(Pt_tmp,Gamma,qf,ds,Sign,0.0,i_panel,j_panel,i,j);
            Pt_tmp=-qc(i,j,1:3); Pt_tmp(1)=-Pt_tmp(1);
            [a1,V4]=wing(Pt_tmp,Gamma,qf,ds,Sign,0.0,i_panel,j_panel,i,j);
        else
            V3(3)=0;V4(3)=0;
        end
        Winduced=V1(3)+V2(3)-V3(3)-V4(3);
        %   add influence of mirror image (ground).
        alpha_induced=-Winduced/Vinf;
        dd(i,j)=rho*dym*Vinf*gamaij*alpha_induced;
        %
        dp(i,j)=dl(i,j)/(ds(i,j,4)*Dynamic_Pressure);
        dly(j)=dly(j)+dl(i,j);
        Lift=Lift+dl(i,j);
        Drag=Drag+dd(i,j);
        Pitching=Pitching+dl(i,j)*(qf(i,j,1)-x(1));
    end
end
Cl=Lift/(Dynamic_Pressure*S);
Cd=Drag/(Dynamic_Pressure*S);
Cm=Pitching/(Dynamic_Pressure*S*c);

%
%        output
%
fprintf('Lift                       = %12.8f \n',Lift);
fprintf('Drag                       = %12.8f \n',Drag);
fprintf('Cl                         = %12.8f \n',Cl);
fprintf('Cd                         = %12.8f \n',Cd);
fprintf('Cm                         = %12.8f \n',Cm);
fprintf('f=Cl/Cd                    = %12.8f \n',Cl/Cd);

%write(6,110)
%line='j  Lift        dCp(1:4)                                             \Gamma_1   \Gamma_1 (2:4)';
%fprintf('%s \n',line);
format='%d  %8.4e   %8.4e   %8.4e   %8.4e   %8.4e   %8.4e   %8.4e   %8.4e   %8.4e \n';
for  j=1:j_panel
    gamma1j(2:i_panel)=diff(Gamma(:,j));
    dlyj(j)=dly(j)/b*j_panel;
%    fprintf(format,j,dlyj,dp(1:4,j),Gamma(1,j),gamma1j(2:4));
end
y=linspace(0,b,j_panel);
figure();
plot(y,dlyj);
title('dLift(y/b)');

figure();
plot(y,dp(1:4,:));
title('dp');

figure();
plot(y,Gamma(1,:));
title('\Gamma');
end

