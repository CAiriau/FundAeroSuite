function table_omega
% Table of the Prandtl-Meyer w.r.t.  Mach number
% C Airiau, avril 2012
k=5
M1=1;
dM=0.01;
n=300;
M2=M1+dM*(n-1);
Mach=M1:dM:M2;
n2=size(Mach);n=n2(2);
if mod(n,k) ~= 0
    disp('bad number in table_omega.m')
end    
omega=Prandtl_Meyer(Mach);
fprintf('Mach dimension ; %i \n',n);
for i=1:k:n
    for j=0:k-2
    fprintf(' %2.2f => %2.2f ',Mach(i+j),omega(i+j));
    end
    fprintf(' %2.2f => %2.2f \n ',Mach(i+k-1),omega(i+k-1));
end
end
