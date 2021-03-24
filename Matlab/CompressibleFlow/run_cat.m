function run_cat(choice)
% main function to run the compressible flow tool CAT
% translation in english : 2020
addpath("cat");

fprintf('number of arguments : %i \n',nargin);
if nargin == 1
	option = choice;
else
	option=0;
end

compressible_flow(option);
display('Normal end of execution of CAT in Matlab');