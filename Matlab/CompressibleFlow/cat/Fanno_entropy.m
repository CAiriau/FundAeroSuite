function Fanno_entropy(Mach,opt)
% Entropy n a Fanno's flow
% opt= 1 : Fanno function, opt= 2 : its derivative wrt Mach

global gam
global r

omega=1.+(gam-1.)/2.*Mach^2;
switch opt
case 1
	tmp=r*(gam+1)/(2.*(gam-1))*log((gam+1)/2.);
    F=r*(log(Mach)-(gam+1)*log(omega)/(2*(gam-1)))+tmp;
otherwise
    F==r*(1.d0-Mach^2)/(Mach*omega);
end 

end

