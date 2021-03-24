function plot_Cp(sc,Cp)
% plot Pressure Coefficient with legend

global name
figure();
plot(sc,Cp,'b-o');
title(name);
xlabel('x/c');
ylabel('Cp');

end
