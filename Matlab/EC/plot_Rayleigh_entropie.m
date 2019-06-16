function plot_Rayleigh_entropie
% dessiner la courbe de l'entropie
 % (s-s_c)/Cp function of Ma
global gam

n=981;
n=201;
Ma_init=0.1;
Ma_end=5;
alpha=(gam+1.0)/2.0;
beta=gam/(gam-1.d0);
 

Ma=linspace(Ma_init,Ma_end,n);
omega=1.0+(gam-1.0)/2.0*Ma.^2;
r_P=(gam+1.0)./(1.0+gam*Ma.^2);
r_T=Ma.^2.*r_P.^2;
r_rho=1.0./(Ma.^2.*r_P);
r_u=1.0./r_rho;
r_Ti=r_T.*omega/alpha;
r_Pi=r_P.*(omega/alpha).^beta;
tmp=2.0*log(Ma)+alpha*log(r_P);

fileID = fopen('Rayleigh.dat','w');
fprintf(fileID,'# \t  M \t (s-sc)/Cp \t\t P/Pc \t T/TC \t\t u/uc \t\t Ti/Tic \t Pi/Pic \t rho/rhoc \n')
for i=1:n    
    fprintf(fileID,'%8.5f \t %8.5f \t %8.5f \t %8.5f \t %8.5f \t %8.5f \t %8.5f \t %8.5f \t  \n',Ma,tmp,r_P,r_T,r_u,r_Ti,r_Pi,r_rho);
end 
fclose(fileID)

figure(1)
title('(s-s_c)/Cp')
plot(Ma,tmp,'b-','LineWidth',2)
xlabel('Mach')
ylabel('(s-s_c)/Cp')
end  
