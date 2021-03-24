function [s,Y,Cp]=grid_generator(vargin) 
% program no. 1: grid generator for 2-d airfoils
%       ----------------------------------------------
%
% Airfoil with van de Vooren and de Jong transformation
%
% Usage :  
%          * grid_generator;
%          * grid_generator(true);
%          * grid_generator([0.1,1.8,10,180]);
%          * grid_generator([e,k,incidence,m]);
%
% Args:    
%          e         : thickness coefficient (lower than 1; ex : 0.1)
%          k         : Trailing edge angle parameter ( < 2, order 2 , ex : 1.8)
%          incidence : angle of attack in °
%          m         : step of theta in degres : 360 / m (must be a multiple of 4!)
%          show_plot : to plot the airfoil and the Cp
%
% 
%
%   this program is an automated complex airfoil transformation of the 
%   type presented by van de vooren and de jong (1970).  the resulting 
%   airfoil may have a non-zero trailing edge angle.  this formulation 
%   is for non-cambered airfoils only (programmed by steven yon, 1989).
%

%%*******************************************************
% Original coding obtained from J. Katz, in fortran 77
% and found in :
% Low-Speed Aerodynamics,Joseph Katz et Allen Plotkin
% Second Edition,Cambridge Aerospace Series, 2001
% 
% translation in Matlab by Christophe Airiau, 2015
% free use of the matlab code, under J. Katz agreement.
%%*******************************************************
tl=1;
profile=1;
switch nargin
    case 0
        disp('ready to start van de vooren transformation');
        e=input('enter thickness coeff. e : ');
        k=input('enter t.e. angle coeff. k : ');
        incidence=input('enter the angle of attack in degrees : ');
        alpha=incidence/180*pi;
        disp('enter number of airfoil panels, m'); 
        disp('with which to model the airfoil');
        disp('(note that m should be an even factor of 360)');
        m=input('enter m (a multiple of 4) =  ');
    case 1
        switch length(vargin)
        case 1 
            e=0.02;
            k=1.85; incidence=0; m=32;show_plot=1;   
            % thickness= 15% e=0.02, k=1.85
        case 5     
        e=vargin(1); k=vargin(2); incidence=vargin(3); m=vargin(4);
        show_plot=vargin(5);
        case 6
        e=vargin(1); k=vargin(2); incidence=vargin(3); m=vargin(4);
        show_plot=vargin(5);
        %profile='parabolic';
        profile=vargin(6);
    otherwise
        disp('bad option');return
    end
end
fid8=fopen('airfoil.dat','w');
fid10=fopen('cp.dat','w');

alpha=incidence/180*pi;
fprintf('e = %f \n',e)


switch profile

case 1
    dtheta=360/m;
    a=2*tl*(e+1)^(k-1)/(2^k); 
    fprintf(' a= %d \n',a);
    %    the do loop will run through the circle plane with
    %    the specified angular interval and transform each 
    %    point to the airfoil plane

    npt=0;
    for i=0:dtheta:360;
        npt=npt+1;
        %fprintf('i = %d \n',i);
        
        if (i <= 0) || (i == 360) 
            x=1; y=0; cp=1;
            fprintf(fid8, '%15.8e  %15.8e ', x,y);
            s(npt)=x;Y(npt)=y;
            if (k == 2) && ((i == 0) || (i == 360) )
                disp('nothing');
            else
                fprintf(fid8, '%15.8e  %15.8e \n ', x,cp);
                Cp(npt)=cp;
            end
        else

            theta=i/180*pi;
            cost=cos(theta);sint=sin(theta);
            r1=abs(a)*sqrt((cost-1)^2+sint^2);
            r2=abs(a)*sqrt((cost-e)^2+sint^2);
            if theta==0 
                theta_1=pi/2;
            else
                theta_1=(atan((a*sint)/(a*(cost-1))))+pi;
            end 

            if(cost-e < 0) && (sint > 0) 
                theta_2=(atan((a*sint)/(a*(cost-e))))+pi;
            elseif (cost-e<0) && (sint<0) 
                theta_2=(atan((a*sint)/(a*(cost-e))))+pi;
            elseif (cost-e>0) && (sint<0) 
                theta_2=(atan((a*sint)/(a*(cost-e))))+2*pi;
            else
                theta_2=(atan((a*sint)/(a*(cost-e))));
            end 

        %   this part computes the transformed positions

            com1=((r1^k)/(r2^(k-1)))/((cos((k-1)*theta_2))^2+(sin((k-1)*theta_2))^2);
            x=com1*(cos(k*theta_1)*cos((k-1)*theta_2) +sin(k*theta_1)*sin((k-1)*theta_2))+tl;
            y=com1*(sin(k*theta_1)*cos((k-1)*theta_2) -cos(k*theta_1)*sin((k-1)*theta_2));
            fprintf(fid8, '%15.8e  %15.8e \n', x,y);
            s(npt)=x;Y(npt)=y;


        %   this part computes the transformed pressure distribution

            a1=cos((k-1)*theta_1)*cos(k*theta_2)+sin((k-1)*theta_1)*sin(k*theta_2);
            b1=sin((k-1)*theta_1)*cos(k*theta_2)-cos((k-1)*theta_1)*sin(k*theta_2);
            c1=(cos(k*theta_2))^2+(sin(k*theta_2))^2;
            p=a*(1-k+k*e);
            d1=a1*(a*cos(theta)-p)-b1*a*sin(theta);
            d2=a1*a*sin(theta)+b1*(a*cos(theta)-p);

            temp=2*c1*(sin(alpha)-sin(alpha-theta))/(d1^2+d2^2);
            com2=temp*(r2^k)/(r1^(k-1));

            vx=d1*sin(theta)+d2*cos(theta);
            vy=-(d1*cos(theta)-d2*sin(theta));

            cp=1-com2^2*(vx^2+vy^2);
            fprintf(fid10, '%15.8e  %15.8e \n', x,cp);Cp(npt)=cp;

        end 
    end
    fprintf('maximun thickness/2 = %f \n',max(abs(Y))) 
    fprintf('min Cp              = %f \n',min(Cp))
    fprintf('alpha               = %f\n',alpha/pi*180)
case 2
    theta=linspace(0,pi,m);
    s=cos(theta)/2;
    Y=e*(1-2*s).*(1+2*s);
    dY=8*e.*s-alpha;
    Cp=-2.*sin(theta).*(4*e+alpha./(1+cos(theta)));
    s=s+0.5;
    for i=1:length(theta)
        fprintf(fid8, '%15.8e  %15.8e \n', s(i),Y(i));
        fprintf(fid10, '%15.8e  %15.8e \n', s(i),Cp(i));
    end
end

fclose(fid8);
fclose(fid10);
if show_plot == 1
    figure(1)
    plot(s,Y,'r-o');
    axis  equal;
    title('airfoil geometry');xlabel('x/chord');
    figure(2)
    plot(s,Cp,'r-o');
    title('Pressure coefficient Cp of reference');xlabel('x/chord');
end

disp('End grid generator')
end

