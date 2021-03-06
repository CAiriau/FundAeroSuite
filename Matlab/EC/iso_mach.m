function iso_mach(choix,option)
%  dessin des iso mach pour des chocs obliques
% C Airiau, avril 2012
if choix==1
    Mach_ini=1.05;
    Mach_fin=8.05;
    dMach=0.1;
    Mach=Mach_ini:dMach:Mach_fin;
else
    Mach=[1.01 1.05 1.1 1.2 1.3  1.4 1.6 1.8 2 2.2 2.5 3 3.5 5 8 10 20  30];
end
figure(1);
fprintf('option = %i \n',option);
if option == 1
    % visualisation en 2D
    for i=1:length(Mach)
    
        fprintf('Mach = %f \n',Mach(i));
        [theta,Sigma]=for_iso_sigma(Mach(i));
        plot(theta,Sigma);
        [x,j]=max(theta);
        y=Sigma(j);
   
        j=50;x=theta(j);y=Sigma(j);
        text(x,y,num2str(Mach(i)));
        xlabel('\theta');ylabel('\sigma');
        hold on
    end
elseif option == 2
    % Visualisation en 3D ...
    for i=1:length(Mach)
         sprintf('Mach = %f \n',Mach(i));
        [theta,Sigma]=for_iso_sigma(Mach(i));
        Z(1:length(Sigma))=1/Mach(i);
        grid on;
        plot3(theta,Sigma,Z);
         xlabel('\theta');ylabel('\sigma');zlabel('1/M_a_v_a_l');
        hold on
    end
    
elseif option == 3
         
    
    [theta,Sigma] = courbe_Mach_aval_unitaire(Mach);
    plot(theta,Sigma,'k-','LineWidth',2);
    grid on;
    xlabel('\theta');ylabel('\sigma');
    [theta,Sigma] = courbe_theta_max(Mach);
    plot(theta,Sigma,'r-','LineWidth',2);
    text(35,20,'red   : max(\theta)');
    text(35,10,'black : M_2=1');
    hold on

end
title('Courbes de choc new');
axis([0 50 0 90]);
grid on;
grid minor;
end
