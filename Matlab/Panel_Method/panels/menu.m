function option=menu
% Menu do define options
%    
%

line='*************************************************';
fprintf('%s\n   Panel method menu \n%s\n',line,line);
fprintf('[1] : grid and airfoil generator \n');
fprintf('[2] : constant strength source   \n');
fprintf('[3] : constant strength doublet  \n');
fprintf('[4] : constant strength vortex   \n');
fprintf('[5] : linear   strength source   \n');
fprintf('[6] : linear   strength vortex   \n');
fprintf('[7] : constant strength doublet potential  \n');
fprintf('[8] : constant strength source doublet potential  \n');
fprintf('[9] : linear   strength doublet potential  \n');
fprintf('[10]: quadratic strength doublet potential  \n');
fprintf('[11]: influence coefficients                \n');
fprintf('[12]: rectangular lifting surface           \n');
fprintf('[13]: lumped vortex                         \n');
fprintf('[14]: unsteady rectangular lifting surface  \n');
fprintf('[15]: Discrete Vortex Method                \n');

option=input('enter your choice please ');

end
