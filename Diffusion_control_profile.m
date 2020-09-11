%This function return the diffusion profile of C by assuming the
%interfacial carbon as equilibrium

function [Diffusion]=Diffusion_control_profile(Xneq,Xpeq,X0,DC,Mob,kafang,Rbcc,SN,distance)
% Diffusion=[Xpm Xip DiffLL Rbcc Xneq Xpeq];
% mol%, mol%, um, um, mol%, mol%

% Xneq=xC_F_eq(i);
% Xpeq=xC_A_eq(i);
% X0=Comp_m(1);
% DC=D_C(i)*1e-12;
% Rbcc=N_pD(l,7);
% SN=length(N_PR{l}(:,1));
% distance=N_p(l,27);
% % distance=N_PR{l}(r,3);

% clear all;
% % Fe-1.0 at.%C as an example
% T=1050; % [K]
% X0=1; % [at.%]
% distance=12.5; % [um]
% V0=4/3*pi*distance^3; % [um^3]
% Xneq=0.0658; % [at.%]
% Xpeq=2.0619; % [at.%]
% DC=1.14e-12; % [m2/s]
% kafang=110; % [J/at.%]
% Mob=5.4e-8; % [mol.m/(J.s)]
% feq=0.1383;
% feq=0.5;
% Vbcc=V0*feq; % [um^3]
% % Abcc=4*pi*(3*Vbcc/(4*pi))^(2/3); % [um^2]
% Rbcc=(3*Vbcc/(4*pi))^(1/3); % [um]
% SN=1;

Vbcc=4/3*pi*Rbcc^3/SN; % [um^3]
syms L;
syms Xm;
fL=Vbcc*(X0-Xneq)*1/SN-(Xpeq-X0)*4*pi*(L^3+5*Rbcc*L^2+10*Rbcc^2*L)/30*1/SN; % mass conservation
sol1=double(vpasolve(fL,L,'random',true));

Lcount=0;
for j=1:length(sol1)
    if isreal(sol1(j))==1 && sol1(j)>0
        if Lcount==0
           DiffLL=sol1(j);
        else if Lcount>1 && sol1(j) < DiffLL
            DiffLL=sol1(j);
            end
        end
        Lcount=Lcount+1;
    end
end
if Lcount==0
   DiffLL=abs(distance-Rbcc);
end
if DiffLL<=(distance-Rbcc)
    softflag=0;
    Xpm=X0;
    Xip=Xpeq;
else
    DiffLL=distance-Rbcc;
    fXm=30*Vbcc*(X0-Xneq)/(4*pi)*1/SN-(9*DiffLL^3*Xm-10*DiffLL^3*X0+DiffLL^3*Xpeq+20*Rbcc^2*Xm*DiffLL- ...
    30*Rbcc^2*X0*DiffLL+10*Rbcc^2*Xpeq*DiffLL+25*DiffLL^2*Rbcc*Xm-30*DiffLL^2*Rbcc*X0+5*DiffLL^2*Rbcc*Xpeq)*1/SN; % mass conservation
    sol2=double(vpasolve(fXm,Xm,'random',true));
    sol2=double(vpasolve(fXm,Xm,[X0 Xpeq],'random',true));
    Mcount=0;
    for j=1:length(sol2)
        if isreal(sol2(j))==1 && sol2(j)>0
            if Mcount==0
               Xpm=sol2(j);
            else if Mcount>1 && sol2(j)>X0 && sol2(j)<Xpeq
                Xpm=sol2(j);
                end
            end
            Mcount=Mcount+1;
        end
    end
    if isempty(sol2)
        Xpm=1/SN*(X0-Xneq*(Rbcc/distance)^3)/(1-(Rbcc/distance)^3);
    end
    Xip=Xpeq;
    softflag=1;
end
Diffusion=[Xpm Xip DiffLL Rbcc Xneq Xpeq softflag]';

% Ngrid=101; % number of grid
% deltaC=X0/Ngrid; % reducing step of C content [at.%]
% iter=1;
% stop=0;
%         
% % calculate the C concentration
% for i=1:Ngrid
%    zz(i)=(i-1)*distance/Ngrid;
%    if zz(i)<Rbcc
%        zC(i)=Xneq; % in ferrite [at.%]
%    else if zz(i)>=Rbcc && zz(i)<=Rbcc+DiffLL
%        zC(i)=Xpm+(Xip-Xpm)*(1-(zz(i)-Rbcc)/DiffLL).^2; % in austenite [at.%]
%        else
%        zC(i)=Xpm;
%        end
%    end
% end
% 
% figure('Name','Carbon concentration profile');
% % plot(zz,zC,'k','LineWidth',2);
% plot(zz,zC,'ko-','LineWidth',2);
% xlabel('z(\mum)');
% ylabel('x_{C} (at.%)');
% box on;
% 
% dlmwrite(strcat(num2str(1),'myfile.txt'),[zz' zC'],'delimiter',' ');

end



