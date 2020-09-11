% This function is to calculate the nucleation rate based on classical nucleation theory
% main Ref: S.E.Offerman et al Science 2002; S.E.Offerman et al Scientific
% Report 2016
% Suitable for Fe-C-Mn system
% CNTt: m-3s-1

function [CNTt EnergyB Freq Scailing]=CNT(Temp,f_bcc,NpotT,deltaGV,t,dt,CR)

% public parameters
h=6.626e-34; % Planck constant [m2.Kg/s] or [J.s]
kB=1.381e-23; % Bolzmann constant [m2.Kg/(s2.K)] or [J/K]
QD=3.93e-19; % activation energy of iron self diffusion [J]
Gs=1e7; % % misfit energy of nucleation [J/m3]
% Gs=2e7; % % misfit energy of nucleation [J/m3], literature shows Gs=1e7 J/m3 at 700C
% T=0K, Gs=2.5*10^8 J/m3, strain =0.03, dV=9 %; T=973 K, dV=1.8%; therefore
% dV (%)=-0.0074*T(K)+9; Gs=Gs0/(e0/e(T))^2;
% Gs=2.5e8./(9./(-0.0074*Temp+9)).^2;

Z=0.065; % Zeldovich factor, nearly constant

% these two parameters should be adjusted to match the experimental data
psi=5e-8; % a factor containing the information of shape of critical nucleus and interfacial energies [J3/m6]; 
          % psi_CF=3.3e-3,psi_LEA=2.1e-6,psi_exp1=5e-8,psi_exp2=1~3e-7

%
SC=0.1; % S/C, scaling factor of potential nucleation sites/nucleation scale
tao=10; % varies in 4~18s according to S.E.Offerman et al Scientific Report 2016

ASC=SC*exp(-tao./t);% pre-factor which could be adjusted by varing SC and tao

if deltaGV<Gs
    Nt=0;
else
    Nt=ASC.*Z*NpotT*(1-f_bcc).*(kB.*Temp./h).*exp(-QD./(kB.*Temp)).*exp(-psi./(kB.*Temp.*(deltaGV-Gs).^2));
end
CNTt=Nt;
EnergyB=exp(-psi./(kB.*Temp.*(deltaGV-Gs).^2));
Freq=(kB.*Temp./h).*exp(-QD./(kB.*Temp));
Scailing=ASC.*Z*NpotT*(1-f_bcc);

end







