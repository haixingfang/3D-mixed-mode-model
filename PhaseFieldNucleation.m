% This function is to calculate the nucleation rate based on Pina's
% arbitrary nucleation mode;
% The aim is to compare with phase field model;
% Ref: M.G.Mecozzi et al, MMTA 2008;
% Suitable for Fe-0.1C-0.49Mn (wt.%) system

function [PFNt]=PhaseFieldNucleation(Temp,CR,dt,Nuclei)
%function [sumN]=CNT(T,f_bcc,Npot,deltaGV,wC_bcc,wC_fcc)

Ts=1092; % [K]
epsT=dt*CR; % [K}
deltaT=24; % [K]
if deltaT==0 && Temp<Ts+epsT && Temp>Ts-epsT
    PFNt=Nuclei/dt;
    else if Temp>Ts-epsT || Temp<(Ts-deltaT)+epsT
    PFNt=0;
    else if Temp<=Ts+epsT && Temp>=Ts-deltaT-epsT
         PFNt=Nuclei/(deltaT/CR);
        end
        end
end

end