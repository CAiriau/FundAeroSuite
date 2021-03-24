function save_data(FileName,x,Cp)
% Write Cp in a text file

%roy = fix(clock);
%date_dat = ['date = ', num2str(roy(1)),'-',num2str(roy(2)),'-',num2str(roy(3)),...
%    ' time = ', num2str(roy(4)),':',num2str(roy(5)),':',num2str(roy(6))]; 

n=length(x);
fid = fopen([FileName '.dat'],'w');
%fprintf(fid,'# File generated at  %s  \n',date_dat);
fprintf(fid,' case : %s\n',FileName);
fprintf(fid,'#          x          Cp \n');
for i = 1:n
    fprintf(fid,'%15.12e %15.12e  %i \n',x(i),Cp(i),i);
end
fclose(fid);
return

end 