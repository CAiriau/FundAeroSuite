function table_omega
% sortie de la table des omega en fonction du Mach
% C Airiau, avril 2012
k=5
M1=1;
dM=0.01;
n=300;
M2=M1+dM*(n-1);
Mach=M1:dM:M2;
n2=size(Mach);n=n2(2);
if mod(n,k) ~= 0
    disp('pas bon nombre')
    
end    
omega=valeur_omega(Mach);
fprintf('dimension de Mach ; %i \n',n);
for i=1:k:n
    for j=0:k-2
    fprintf(' %2.2f => %2.2f ',Mach(i+j),omega(i+j));
    end
    fprintf(' %2.2f => %2.2f \n ',Mach(i+k-1),omega(i+k-1));
end
end
