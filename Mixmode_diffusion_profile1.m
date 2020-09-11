%This function return the diffusion profile of C according to C.Bos and
%J.Sietsma, 2007

function [Diffusion]=Mixmode_diffusion_profile1(Xneq,Xpeq,X0,DC,Mob,kafang,Rbcc,SN,distance)
% Diffusion=[Xpm Xip DiffLL Rbcc Xneq Xpeq];
% mol%, mol%, um, um, mol%, mol%

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
% feq=0.3;
% Vbcc=V0*feq; % [um^3]
% % Abcc=4*pi*(3*Vbcc/(4*pi))^(2/3); % [um^2]
% Rbcc=(3*Vbcc/(4*pi))^(1/3); % [um]
% SN=1;

Vbcc=4/3*pi*Rbcc^3/SN; % [um^3]
Z=2*DC/(Mob*kafang)*1e6; % [um*at.%]
FF=@(x)mixmode_fun(x,X0,Xneq,Xpeq,Rbcc,Z,SN);
opts=optimset('Display','iter','Algorithm',{'levenberg-marquardt',0.005},'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-6);% option settings for fsolve
x0(1)=Xpeq; % initial Xip [at.%]
% x0(2)=distance-Rbcc; % initial diffusion length [um]
x0(2)=0.05;
[sol1,fval,exitflag]=fsolve(FF,x0,opts);
Xip=sol1(1);
DiffLL=sol1(2);
Xpm=X0;

if DiffLL<=(distance-Rbcc)
    softflag=0;
else
    SS=@(xx)mixmode_fun_softimpingement(xx,X0,Xneq,Xpeq,Rbcc,Z,distance,SN)
    opts=optimset('Display','iter','Algorithm',{'levenberg-marquardt',0.005},'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-6);% option settings for fsolve
    xx0(1)=(X0-(Rbcc/distance)^3*Xneq)/(1-(Rbcc/distance)^3); % initial Xip [at.%]
    xx0(2)=distance-Rbcc; % initial diffusion length [um]
    xx0(3)=(X0-(Rbcc/distance)^3*Xneq)/(1-(Rbcc/distance)^3); % initial central C content[at.%]
    [sol2,fval,exitflag]=fsolve(SS,xx0,opts);
    Xip=sol2(1);
    DiffLL=sol2(2);
    Xpm=sol2(3);
    softflag=1;
end
Diffusion=[Xpm Xip DiffLL Rbcc Xneq Xpeq softflag]';
% Ngrid=101; % number of grid
% deltaC=X0/Ngrid; % reducing step of C content [at.%]
% iter=1;
% stop=0;
        
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



