
% Function of magnetic configuration for ferrite_3d_model
% This function is to calculate the magnetic configuration from the
% microstruction derived from ferrite_3d_model.m by using fgoalattain
% function 

% Note: for all matrix vectors: colume 1, 2, 3 denotes x,y,z components
% respectively. Unit of magnetic field is either mA/m or T except for special
% cases
% Mk is in [A/m]

function [Mag_parameter0,kesa0,md_cos0,nd_angle0,nd_delta0,aax0]=Magnetic_calculation_modify_try(alpha0,beta0,MA_theta0,MA_phi0,x,xx,numb,f,temperature,Lb,Tc)
% Input:
% alpha0 and beta0 is the assigned initial values for magnetization unit vector describe in polar coordinate
% MA_alpha0 and MA_beta0 is the assigned initial values for anistropy axis unit vector describe in polar coordinate
% x is the N0*12 arrat that contains all the microstructure information
% f is the volume fraction of ferromagnetic phase (ferrite)
% numb is the number of nucleation site including the inactive nucleation site
% temperature is used to calculate the saturation magnetization

% Output:
% Mag_parameter0 contains all the magnetic configuration
% kesa0 is the correlation function
% md_cos0 is the mean direction cosine, a measure for the magnetic texture in xyz direction
% nd_angle0 is the rotation angle of the polarization vector
% nd_delta0 is the mean size derived from ND simulation
% aax0 is the number of magnetic field point

% clear all;
% Npot=30;
% numb=20;
% f=0.15;
% temperature=1030;
% 
% alpha0=2*pi*rand(numb,1);
% beta0=pi*rand(numb,1);
% MA_theta0  = 2*pi*rand(1,Npot);
% MA_phi0    = pi*rand(1,Npot);
% x=xlsread('D:\hfang\Desktop\results.xlsx','sheet1','A2:L21');% Read the matrix of microstructure

% Some basic parameters
dmax=5.2;            % Maximum ratio of distance ij to the radius of particle j reflecting the influence range of dipole field from j
inter_alpha=1/3;     % Inter domain coupling constant representing the mean field interaction
Ns=0.5;              % Shape anistropy factor between -0.5 and 1
u0=4*pi*1e-7;        % Vaccum permeability [H/m]
%Lb=126;              % Lb the side of the cubix box, same as that in ferrite_3d_model [um]
distribution_a=0.5;  % The constant in the Gaussian distribution for the delta function
syms fai;
syms fai_c;

% Following is to calculate the saturate magnetization Bs at a certain temperature
% Reference: A.S. Arrott and B. Heinrich, J. Appl. Phys. 52 (1981) 2113.
%Tc=1043;   % Curie temperature [K]
B0=2.20;   % magnetic field at 0 K [T]
Arrott_beta = 0.368; %critical exponent beta
Arrott_A = 0.110;
Arrott_C = 0.129;
tau=temperature/Tc; % reduced temperature
if temperature<Tc
   Bs=B0*(1-tau)^Arrott_beta/(1-Arrott_beta*tau+Arrott_A*tau^(3/2)-Arrott_C*tau^(7/2));% Saturation magnetic induction [T]
   else
   Bs=0;
end
Ms=Bs/u0;% Saturation magnetization [A/m]

% Assign the magnetic anisotropy vector for the potential particle at time t
MA_theta  = MA_theta0(1,1:numb);
MA_phi    = MA_phi0(1,1:numb);
Mag_anist(:,1) = sin(MA_phi).*cos(MA_theta);
Mag_anist(:,2) = sin(MA_phi).*sin(MA_theta);
Mag_anist(:,3) = cos(MA_phi);

% Assign the initial value for magnetization of each ferrite grain
alpha_m=alpha0(1:numb,1);
beta_m=beta0(1:numb,1);

stop_flag=2;% Sign for the stop of the calculation,2 means positive calculation,1 means negative calculation,0 means stop
direction=1;% Assign the initial calculation direction: from positive to negative
while stop_flag~=0
% Assign the applied field vector
aax0=1;
for aa=1:aax0
Ba(aa,1)=0;
Ba(aa,2)=0;
if aax0>1
  if direction==1            % Calculate from positive to negative it is 1,otherwise -1
     Ba(aa,3)=0.62-aa*1.2/(aax0-1);
  else
     Ba(aa,3)=-0.62+aa*1.2/(aax0-1);% Calculate from negative to positive it is -1,otherwise 1
  end
else
   Ba(aa,3)=0.3;
end
end
Ha=Ba./u0; % Applied field in [A/m]

for aa=1:aax0
  
% Following is to calculate the magnetic field direction for each particle
%Assign the initial value of total average magnetization

stop_cycle=0;% Stop sign for the main loop
cc=1;% Restore the number of the main loop
while stop_cycle~=1;
      sum(aa)=0;
      M(1)=0;M(2)=0;M(3)=0; % Assign the magnetization for each particle and calculate the total average magnetizatino first
      inactive=0;
      for k=1:numb
          if cc==1&&aa==1
             Mk(k,1)=Ms*cos(alpha_m(k))*sin(beta_m(k));
             Mk(k,2)=Ms*sin(alpha_m(k))*sin(beta_m(k));
             Mk(k,3)=Ms*cos(beta_m(k));
          else if cc==1&&aa>1
             Mk(k,1)=Magnetization(aa-1,k,1);
             Mk(k,2)=Magnetization(aa-1,k,2);
             Mk(k,3)=Magnetization(aa-1,k,3);
             else if cc>1
             Mk(k,1)=Mkk(k,1);
             Mk(k,2)=Mkk(k,2);
             Mk(k,3)=Mkk(k,3);
                  end
              end
          end
         % The magnetization of the inative nucleated site should be zero
         if x(k,5)==0
            Mk(k,:)=0;
            inactive=inactive+1; % Count of the number of inactive nucleation site
         end
          
         M(1)=M(1)+Mk(k,1)*x(k,10)/(Lb^3);% Overall mean magnetization in X direction
         M(2)=M(2)+Mk(k,2)*x(k,10)/(Lb^3);% Overall mean magnetization in Y direction
         M(3)=M(3)+Mk(k,3)*x(k,10)/(Lb^3);% Overall mean magnetization in Z direction
      end
      
for i=1:numb
    Hdd(aa,i,1)=0;% Assign the initial value of dipole field Hdd for particle i
    Hdd(aa,i,2)=0;
    Hdd(aa,i,3)=0;
    for j=1:numb
        if i~=j
           d(j,1)=x(i,2)-x(j,2);
           d(j,2)=x(i,3)-x(j,3);
           d(j,3)=x(i,4)-x(j,4); % Distance vector of particle i and j
           d_MXY(j,1)=x(i,2)-xx(j+numb,1);        
           d_MXY(j,2)=x(i,3)-xx(j+numb,2);
           d_MXY(j,3)=x(i,4)-xx(j+numb,3); % Distance vector of particle i and j mirror point XY
           d_MXZ(j,1)=x(i,2)-xx(j+numb*2,1);
           d_MXZ(j,2)=x(i,3)-xx(j+numb*2,2);
           d_MXZ(j,3)=x(i,4)-xx(j+numb*2,3); % Distance vector of particle i and j mirror point XZ
           d_MYZ(j,1)=x(i,2)-xx(j+numb*3,1);
           d_MYZ(j,2)=x(i,3)-xx(j+numb*3,2);
           d_MYZ(j,3)=x(i,4)-xx(j+numb*3,3); % Distance vector of particle i and j mirror point YZ           
           d_MO(j,1)=x(i,2)-xx(j+numb*4,1);
           d_MO(j,2)=x(i,3)-xx(j+numb*4,2);
           d_MO(j,3)=x(i,4)-xx(j+numb*4,3); % Distance vector of particle i and j mirror point O           

           dis(i,j)=sqrt(d(j,1)^2+d(j,2)^2+d(j,3)^2); % Distance between grain 1 and 2
           dis_MXY(i,j)=sqrt(d_MXY(j,1)^2+d_MXY(j,2)^2+d_MXY(j,3)^2); % Distance between grain 1 and 2 mirrored by XY
           dis_MXZ(i,j)=sqrt(d_MXZ(j,1)^2+d_MXZ(j,2)^2+d_MXZ(j,3)^2); % Distance between grain 1 and 2 mirrored by XZ
           dis_MYZ(i,j)=sqrt(d_MYZ(j,1)^2+d_MYZ(j,2)^2+d_MYZ(j,3)^2); % Distance between grain 1 and 2 mirrored by YZ
           dis_MO(i,j)=sqrt(d_MO(j,1)^2+d_MO(j,2)^2+d_MO(j,3)^2); % Distance between grain 1 and 2 mirrored by O
           
           dis_min(i,j)=min([dis(i,j),dis_MXY(i,j),dis_MXZ(i,j),dis_MYZ(i,j),dis_MO(i,j)]); % Get the minimum dista
           pointer=find([dis(i,j),dis_MXY(i,j),dis_MXZ(i,j),dis_MYZ(i,j),dis_MO(i,j)]==dis_min(i,j));
           
           ratio(i,j)=(dis_min(i,j)-x(j,8))/x(j,8);
           
               p(j,1)=x(j,10)*Mk(j,1)*1e-18;
               p(j,2)=x(j,10)*Mk(j,2)*1e-18;
               p(j,3)=x(j,10)*Mk(j,3)*1e-18;% Magnetic dipole moment vector of particle j [Am^2]
                    
           if dis_min(i,j)/x(j,8)<=dmax    % The dipole influence is limited
             if pointer==1
               Hd(j,1)=1/(4*pi)*((3*d(j,1)*(d(j,1)*p(j,1)+d(j,2)*p(j,2)+d(j,3)*p(j,3))-dis(i,j)^2*p(j,1))/(dis(i,j)^5*1e-18));
               Hd(j,2)=1/(4*pi)*((3*d(j,2)*(d(j,1)*p(j,1)+d(j,2)*p(j,2)+d(j,3)*p(j,3))-dis(i,j)^2*p(j,2))/(dis(i,j)^5*1e-18));
               Hd(j,3)=1/(4*pi)*((3*d(j,3)*(d(j,1)*p(j,1)+d(j,2)*p(j,2)+d(j,3)*p(j,3))-dis(i,j)^2*p(j,3))/(dis(i,j)^5*1e-18));% [A/m]
             else if pointer==2
               Hd(j,1)=1/(4*pi)*((3*d_MXY(j,1)*(d_MXY(j,1)*p(j,1)+d_MXY(j,2)*p(j,2)+d_MXY(j,3)*p(j,3))-dis_MXY(i,j)^2*p(j,1))/(dis_MXY(i,j)^5*1e-18));
               Hd(j,2)=1/(4*pi)*((3*d_MXY(j,2)*(d_MXY(j,1)*p(j,1)+d_MXY(j,2)*p(j,2)+d_MXY(j,3)*p(j,3))-dis_MXY(i,j)^2*p(j,2))/(dis_MXY(i,j)^5*1e-18));
               Hd(j,3)=1/(4*pi)*((3*d_MXY(j,3)*(d_MXY(j,1)*p(j,1)+d_MXY(j,2)*p(j,2)+d_MXY(j,3)*p(j,3))-dis_MXY(i,j)^2*p(j,3))/(dis_MXY(i,j)^5*1e-18));% [A/m]               
               
             else if pointer==3
               Hd(j,1)=1/(4*pi)*((3*d_MXZ(j,1)*(d_MXZ(j,1)*p(j,1)+d_MXZ(j,2)*p(j,2)+d_MXZ(j,3)*p(j,3))-dis_MXZ(i,j)^2*p(j,1))/(dis_MXZ(i,j)^5*1e-18));
               Hd(j,2)=1/(4*pi)*((3*d_MXZ(j,2)*(d_MXZ(j,1)*p(j,1)+d_MXZ(j,2)*p(j,2)+d_MXZ(j,3)*p(j,3))-dis_MXZ(i,j)^2*p(j,2))/(dis_MXZ(i,j)^5*1e-18));
               Hd(j,3)=1/(4*pi)*((3*d_MXZ(j,3)*(d_MXZ(j,1)*p(j,1)+d_MXZ(j,2)*p(j,2)+d_MXZ(j,3)*p(j,3))-dis_MXZ(i,j)^2*p(j,3))/(dis_MXZ(i,j)^5*1e-18));% [A/m]               
               
             else if pointer==4
               Hd(j,1)=1/(4*pi)*((3*d_MYZ(j,1)*(d_MYZ(j,1)*p(j,1)+d_MYZ(j,2)*p(j,2)+d_MYZ(j,3)*p(j,3))-dis_MYZ(i,j)^2*p(j,1))/(dis_MYZ(i,j)^5*1e-18));
               Hd(j,2)=1/(4*pi)*((3*d_MYZ(j,2)*(d_MYZ(j,1)*p(j,1)+d_MYZ(j,2)*p(j,2)+d_MYZ(j,3)*p(j,3))-dis_MYZ(i,j)^2*p(j,2))/(dis_MYZ(i,j)^5*1e-18));
               Hd(j,3)=1/(4*pi)*((3*d_MYZ(j,3)*(d_MYZ(j,1)*p(j,1)+d_MYZ(j,2)*p(j,2)+d_MYZ(j,3)*p(j,3))-dis_MYZ(i,j)^2*p(j,3))/(dis_MYZ(i,j)^5*1e-18));% [A/m]               

             else if pointer==5
               Hd(j,1)=1/(4*pi)*((3*d_MO(j,1)*(d_MO(j,1)*p(j,1)+d_MO(j,2)*p(j,2)+d_MO(j,3)*p(j,3))-dis_MO(i,j)^2*p(j,1))/(dis_MO(i,j)^5*1e-18));
               Hd(j,2)=1/(4*pi)*((3*d_MO(j,2)*(d_MO(j,1)*p(j,1)+d_MO(j,2)*p(j,2)+d_MO(j,3)*p(j,3))-dis_MO(i,j)^2*p(j,2))/(dis_MO(i,j)^5*1e-18));
               Hd(j,3)=1/(4*pi)*((3*d_MO(j,3)*(d_MO(j,1)*p(j,1)+d_MO(j,2)*p(j,2)+d_MO(j,3)*p(j,3))-dis_MO(i,j)^2*p(j,3))/(dis_MO(i,j)^5*1e-18));% [A/m]               
                 end
                 end
                 end
                 end
             end
               Hdd(aa,i,1)=Hdd(aa,i,1)+Hd(j,1);
               Hdd(aa,i,2)=Hdd(aa,i,2)+Hd(j,2);
               Hdd(aa,i,3)=Hdd(aa,i,3)+Hd(j,3);% Local field vector plus the contribution from dipole field
           end
        end
    end
    % Interactions including dipole-dipole field and mean field
    Hk(aa,i,1)=Ha(aa,1)+Hdd(aa,i,1)+inter_alpha*M(1);% Ha+Hdd+Hf
    Hk(aa,i,2)=Ha(aa,2)+Hdd(aa,i,2)+inter_alpha*M(2);% Ha+Hdd+Hf
    Hk(aa,i,3)=Ha(aa,3)+Hdd(aa,i,3)+inter_alpha*M(3);% Ha+Hdd+Hf
    
    % noninteractions
%     Hk(aa,i,1)=Ha(aa,1);% Local field equals to external applied field
%     Hk(aa,i,2)=Ha(aa,2);% Local field equals to external applied field
%     Hk(aa,i,3)=Ha(aa,3);% Local field equals to external applied field
      
    % The local field of the inative nucleated site should be zero
    if x(i,5)==0
       Hk(aa,i,:)=0.0001;
    end    
    Hkk(aa,i)=sqrt(Hk(aa,i,1)^2+Hk(aa,i,2)^2+Hk(aa,i,3)^2); % Length of the vector Hk
    theta(aa,i)=acos((Mag_anist(i,1)*Hk(aa,i,1)+Mag_anist(i,2)*Hk(aa,i,2)+Mag_anist(i,3)*Hk(aa,i,3))/(1*Hkk(aa,i))); % The angle between easy axis and Hk
    if aa==1
    sw(i)=0;
    end
    if theta(aa,i)>pi/2
        theta(aa,i)=pi-theta(aa,i);% Always keep this angle between 0 and pi/2, meaning the easy axis direction change with the field direction
        Mag_anist(i,1)=-Mag_anist(i,1);
        Mag_anist(i,2)=-Mag_anist(i,2);
        Mag_anist(i,3)=-Mag_anist(i,3);
        sw(i)=sw(i)+1; % To trace the times of magnetization swithing 
    end
    theta_t(aa,i)=nthroot(tan(theta(aa,i)),3);
    K_const=0.5*u0*Ns*Ms^2; % Anistropy constant [J/m3]
    H_Am=K_const/(u0*Ms); % Imagined anistropy field [A/m]
    hh(aa,i)=u0*Ms*Hkk(aa,i)/(2*K_const);   
%        Hk_initial(i,1)=0;
%        Hk_initial(i,2)=0;
%        Hk_initial(i,3)=1;
       Hk_initial(i,1) = sin(MA_phi(i))*cos(MA_theta(i));
       Hk_initial(i,2) = sin(MA_phi(i))*sin(MA_theta(i));
       Hk_initial(i,3) = abs(cos(MA_phi(i)));    
       Hkk_initial(i)=1;    
    rotation_angle(aa,i)=acos((Hk(aa,i,1)*Hk_initial(i,1)+Hk(aa,i,2)*Hk_initial(i,2)+Hk(aa,i,3)*Hk_initial(i,3))/(Hkk(aa,i)*Hkk_initial(i)));
    if rotation_angle(aa,i)>=pi/2
        hh(aa,i)=-hh(aa,i);
    end
    hs(aa,i)=sqrt(1-theta_t(aa,i)^2+theta_t(aa,i)^4)/(1+theta_t(aa,i)^2);  
    
    sol_flag=0; % The flag for evaluating if the second derivative >0 is meet or not
    uu=1;
    while sol_flag~=1
    
    if direction==1
       if abs(hh(aa,i))>hs(aa,i)&&hh(aa,i)>=0
          sol1=vpasolve((1/2*sin(2*(fai-theta(aa,i)))+(hh(aa,i))*sin(fai)),fai,[0 pi],'random',true);
          fai_solution(aa,i)=double(sol1);
          E_total(aa,i)=1/4-1/4*cos(2*(theta(aa,i)-fai_solution(aa,i)))-hh(aa,i)*cos(fai_solution(aa,i));% total energy
          c(aa,i)=cos(2*(fai_solution(aa,i)-theta(aa,i)))+(hh(aa,i))*cos(fai_solution(aa,i));% second-derivative should be > 0
          if c(aa,i)>0
             sol_flag=1;
          end
       else if abs(hh(aa,i))<=hs(aa,i)&&hh(aa,i)>=0
          sol1=vpasolve((1/2*sin(2*(fai-theta(aa,i)))+(hh(aa,i))*sin(fai)),fai,[0 pi],'random',true); 
          fai_solution(aa,i)=double(sol1);
          E_total(aa,i)=1/4-1/4*cos(2*(theta(aa,i)-fai_solution(aa,i)))-hh(aa,i)*cos(fai_solution(aa,i));% total energy
          c(aa,i)=cos(2*(fai_solution(aa,i)-theta(aa,i)))+(hh(aa,i))*cos(fai_solution(aa,i));% second-derivative should be > 0
          if c(aa,i)>0
             sol_flag=1;
          end
       else if abs(hh(aa,i))<=hs(aa,i)&&hh(aa,i)<0
          sol1=vpasolve((1/2*sin(2*(fai-theta(aa,i)))+(hh(aa,i))*sin(fai)),fai,[0 pi],'random',true);
          fai_solution(aa,i)=double(sol1);
          E_total(aa,i)=1/4-1/4*cos(2*(theta(aa,i)-fai_solution(aa,i)))-hh(aa,i)*cos(fai_solution(aa,i));% total energy
          c(aa,i)=cos(2*(fai_solution(aa,i)-theta(aa,i)))+(hh(aa,i))*cos(fai_solution(aa,i));% second-derivative should be > 0
          if c(aa,i)>0
             sol_flag=1;
          end
               else if abs(hh(aa,i))>hs(aa,i)&&hh(aa,i)<0
                  sol1=vpasolve((1/2*sin(2*(fai-theta(aa,i)))+(hh(aa,i))*sin(fai)),fai,[pi 2*pi],'random',true); 
                  fai_solution(aa,i)=double(sol1);
                  E_total(aa,i)=1/4-1/4*cos(2*(theta(aa,i)-fai_solution(aa,i)))-hh(aa,i)*cos(fai_solution(aa,i));% total energy
                  c(aa,i)=cos(2*(fai_solution(aa,i)-theta(aa,i)))+(hh(aa,i))*cos(fai_solution(aa,i));% second-derivative should be > 0
                  if c(aa,i)>0
                     sol_flag=1;
                  end
                   end
              end
           end
       end
    end
       
     if direction==-1
       if abs(hh(aa,i))>hs(aa,i)&&hh(aa,i)<0
          sol1=vpasolve((1/2*sin(2*(fai-theta(aa,i)))+(hh(aa,i))*sin(fai)),fai,[pi 2*pi],'random',true);
          fai_solution(aa,i)=double(sol1);
          c(aa,i)=cos(2*(fai_solution(aa,i)-theta(aa,i)))+(hh(aa,i))*cos(fai_solution(aa,i));% second-derivative should be > 0
          if c(aa,i)>0
             sol_flag=1;
          end
       else if abs(hh(aa,i))<=hs(aa,i)&&hh(aa,i)<0
          sol1=vpasolve((1/2*sin(2*(fai-theta(aa,i)))+(hh(aa,i))*sin(fai)),fai,[pi 2*pi],'random',true); 
          fai_solution(aa,i)=double(sol1);
          c(aa,i)=cos(2*(fai_solution(aa,i)-theta(aa,i)))+(hh(aa,i))*cos(fai_solution(aa,i));% second-derivative should be > 0
          if c(aa,i)>0
             sol_flag=1;
          end
          else if abs(hh(aa,i))<=hs(aa,i)&&hh(aa,i)>=0
            sol1=vpasolve((1/2*sin(2*(fai-theta(aa,i)))+(hh(aa,i))*sin(fai)),fai,[pi 2*pi],'random',true);
            fai_solution(aa,i)=double(sol1);
            c(aa,i)=cos(2*(fai_solution(aa,i)-theta(aa,i)))+(hh(aa,i))*cos(fai_solution(aa,i));% second-derivative should be > 0
            if c(aa,i)>0
                sol_flag=1;
            end
               else if abs(hh(aa,i))>hs(aa,i)&&hh(aa,i)>=0
                      sol1=vpasolve((1/2*sin(2*(fai-theta(aa,i)))+(hh(aa,i))*sin(fai)),fai,[0 pi],'random',true);
                      fai_solution(aa,i)=double(sol1);
                      c(aa,i)=cos(2*(fai_solution(aa,i)-theta(aa,i)))+(hh(aa,i))*cos(fai_solution(aa,i));% second-derivative should be > 0
                         if c(aa,i)>0
                            sol_flag=1;
                         end
                   end
              end
           end
       end
     end
     uu=uu+1;
    end
         
       Eh(aa,i)=-u0*Ms*cos(fai_solution(aa,i))*Hkk(aa,i); %Local field energy density Eh [J/m3] 
       Ea(aa,i)=K_const*(sin(theta(aa,i)-fai_solution(aa,i)))^2; % Anistropy energy density Ea [J/m3]
       Obj(aa,i)=H_Am*sin(2*(theta(aa,i)-fai_solution(aa,i)))-Hkk(aa,i)*sin(fai_solution(aa,i));% First derivative should be = 0
       sum(aa)=sum(aa)+Eh(aa,i)+Ea(aa,i);   
end

MM(1)=0;MM(2)=0;MM(3)=0;
for i=1:numb
    if direction==1
    if abs(hh(aa,i))>hs(aa,i)&&hh(aa,i)>=0
       fai_angle(aa,i)=fai_solution(aa,i);
       theta_angle(aa,i)=theta(aa,i)-fai_solution(aa,i);
    else if abs(hh(aa,i))<=hs(aa,i)&&hh(aa,i)>=0
       fai_angle(aa,i)=fai_solution(aa,i);
       theta_angle(aa,i)=theta(aa,i)-fai_solution(aa,i);
        else if abs(hh(aa,i))<=hs(aa,i)&&hh(aa,i)<0
                fai_angle(aa,i)=pi-fai_solution(aa,i);
                theta_angle(aa,i)=pi-(theta(aa,i)-fai_solution(aa,i));
            else if abs(hh(aa,i))>hs(aa,i)&&hh(aa,i)<0
                        fai_angle(aa,i)=fai_solution(aa,i)-pi;
                        theta_angle(aa,i)=theta(aa,i)-(fai_solution(aa,i)-pi);
                    end
            end
        end
    end
    else
      if abs(hh(aa,i))>hs(aa,i)&&hh(aa,i)<0
       fai_angle(aa,i)=fai_solution(aa,i)-pi;
       theta_angle(aa,i)=theta(aa,i)-(fai_solution(aa,i)-pi);
    else if abs(hh(aa,i))<=hs(aa,i)&&hh(aa,i)<0
       fai_angle(aa,i)=fai_solution(aa,i)-pi;
       theta_angle(aa,i)=theta(aa,i)-(fai_solution(aa,i)-pi);
        else if abs(hh(aa,i))<=hs(aa,i)&&hh(aa,i)>=0
                fai_angle(aa,i)=2*pi-fai_solution(aa,i);
                theta_angle(aa,i)=fai_solution(aa,i)-theta(aa,i);
            else if abs(hh(aa,i))>hs(aa,i)&&hh(aa,i)>=0
                        fai_angle(aa,i)=fai_solution(aa,i);
                        theta_angle(aa,i)=theta(aa,i)-fai_solution(aa,i);
                    end
            end
        end
            end
    end   
    ff=@(N)magnetization_vector(N,Mag_anist(i,:),Hk(aa,i,:),fai_angle(aa,i),theta_angle(aa,i));
    opts=optimset('Display','iter','Algorithm',{'levenberg-marquardt',0.005},'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-6);% option settings for fsolve
    out0(1)=1*2*(rand(1,1)-0.5);
    out0(2)=1*2*(rand(1,1)-0.5);
    [out,fval3,exitflag3]=fsolve(ff,out0,opts);
    xx3(aa,i,1)=out(1);
    xx3(aa,i,2)=out(2);
    judge3(cc,i)=exitflag3;
      
    Mkk(i,1)=(out(1)*Hk(aa,i,1)/Hkk(aa,i)+out(2)*Mag_anist(i,1))*Ms;
    Mkk(i,2)=(out(1)*Hk(aa,i,2)/Hkk(aa,i)+out(2)*Mag_anist(i,2))*Ms;
    Mkk(i,3)=(out(1)*Hk(aa,i,3)/Hkk(aa,i)+out(2)*Mag_anist(i,3))*Ms;
    
    if Mkk(i,3)>Ms
       Mkk(i,3)=Ms;
    end
    beta(cc,i)=acos(Mkk(i,3)/Ms);
    alpha(cc,i)=atan2(Mkk(i,2)/Ms,Mkk(i,1)/Ms);
    if alpha(cc,i)<0
       alpha(cc,i)=alpha(cc,i)+pi;
    end

    MM(1)=MM(1)+Mkk(i,1)*x(i,10)/(Lb^3);% Updated overall mean magnetization in X direction
    MM(2)=MM(2)+Mkk(i,2)*x(i,10)/(Lb^3);% Updated overall mean magnetization in Y direction
    MM(3)=MM(3)+Mkk(i,3)*x(i,10)/(Lb^3);% Updated overall mean magnetization in Z direction
end
 
  M_total(cc,1)=MM(1);
  M_total(cc,2)=MM(2);
  M_total(cc,3)=MM(3);
  if cc>1
    %flag=abs((sqrt(M_total(cc,1)^2+M_total(cc,2)^2+M_total(cc,3)^2)-sqrt(M_total(cc-1,1)^2+M_total(cc-1,2)^2+M_total(cc-1,3)^2))/sqrt(M_total(cc,1)^2+M_total(cc,2)^2+M_total(cc,3)^2));
    flag1(aa)=(M_total(cc,1)-M_total(cc-1,1))/M_total(cc,1);
    flag2(aa)=(M_total(cc,2)-M_total(cc-1,2))/M_total(cc,2);
    flag3(aa)=(M_total(cc,3)-M_total(cc-1,3))/M_total(cc,3);
    flag4(aa)=abs(sqrt(M_total(cc,1)^2+M_total(cc,2)^2+M_total(cc,3)^2)-sqrt(M_total(cc-1,1)^2+M_total(cc-1,2)^2+M_total(cc-1,3)^2))/sqrt(M_total(cc-1,1)^2+M_total(cc-1,2)^2+M_total(cc-1,3)^2);
    %if flag1<=0.1&&flag2<=0.1&&flag3<=0.1&&flag4<=0.1
    if flag4(aa)<=0.005||cc>20
    %if flag4(aa)<=0.005
       stop_cycle=1;
       for i=1:numb
       Magnetization(aa,i,1)=Mkk(i,1);
       Magnetization(aa,i,2)=Mkk(i,2);
       Magnetization(aa,i,3)=Mkk(i,3);
       end
    end
  end

cc=cc+1
end

% Plot the magnetization vector as arrows for each particle
if aax0>1
if direction==1
    if aa==1||aa==fix(aax0/2)||aa==aax0
        subplot(2,2,3);
        quiver3(x(:,2),x(:,3),x(:,4),Mk(:,1),Mk(:,2),Mk(:,3),'AutoScale','on','LineWidth',1.5,'MaxHeadSize',0.1);
        set(gca,'FontSize',14);
        view(-45,45);
        xlabel('x','FontSize',20);
        ylabel('y','FontSize',20);
        zlabel('z','FontSize',20);
        hold on
        if aa==aax0
           quiver3(x(:,2),x(:,3),x(:,4),Mag_anist(:,1),Mag_anist(:,2),Mag_anist(:,3),'r:','AutoScale','on','LineWidth',1.5,'ShowArrowHead','off'); 
           legend(sprintf('Ba = %f T', Ba(1,3)),sprintf('Ba = %f T', Ba(fix(aax0/2),3)),sprintf('Ba = %f T', Ba(aax0,3)),'Easy axis','FontSize',20);
           title('Magnetization vector','FontSize',20);
        end
    end
end
end
% Characteristics of magnetic structure
    sum_nx=0;sum_ny=0;sum_nz=0;
    sum_mx=0;sum_my=0;sum_mz=0; 
    sum_V=0;
    nnx(aa)=0;nny(aa)=0;nnz(aa)=0;
for i=1:numb    
        n(i,:)=Mkk(i,:)./(sqrt(Mkk(i,1)^2+Mkk(i,2)^2+Mkk(i,3)^2)); 
        
        sum_mx=sum_mx+Mkk(i,1)*x(i,10); 
        sum_my=sum_my+Mkk(i,2)*x(i,10); 
        sum_mz=sum_mz+Mkk(i,3)*x(i,10); 
        
        sum_nx=sum_nx+n(i,1)^2*x(i,10); 
        sum_ny=sum_ny+n(i,2)^2*x(i,10); 
        sum_nz=sum_nz+n(i,3)^2*x(i,10);
        
        nnx(aa)=nnx(aa)+n(i,1)^2;   % summation of the x-direction indices
        nny(aa)=nny(aa)+n(i,2)^2;   % summation of the y-direction indices
        nnz(aa)=nnz(aa)+n(i,3)^2;   % summation of the z-direction indices
        
        sum_V=sum_V+x(i,10); 
end  
nx(aa)=sum_nx/sum_V;
ny(aa)=sum_ny/sum_V;
nz(aa)=sum_nz/sum_V;
mx(aa)=sum_mx/(sum_V*Ms);
my(aa)=sum_my/(sum_V*Ms);
mz(aa)=sum_mz/(sum_V*Ms);
if mz(aa)>0
    mm(aa)=sqrt(mx(aa)^2+my(aa)^2+mz(aa)^2);
else
    mm(aa)=-sqrt(mx(aa)^2+my(aa)^2+mz(aa)^2);
end

% calculate the ND characteristics from the magnetic configuration

% Some basic parameters
gamma=1.83*1e8;   % Gyromagnetic factor [s-1T-1]
wavelength=0.18;  % Usually use the thermal neutrons [nm]
c=2.15*1e29*(wavelength*1e-9)^2; % A constant [T-2m-4]
L=Lb;            % The thickness of the sample [um]
c2=(4*pi*f^2/81)^(1/3);
average4_R=0;
average3_R=0;
    for i=1:numb
        average4_R=average4_R+x(i,8)^4; % Average value of R^4 [um^4]
        average3_R=average3_R+x(i,8)^3; % Average value of R^3 [um^3]
    end
average4_R=average4_R/numb; %[um^4]
average3_R=average3_R/numb; %[um^3]
    
for aa=1:aax0
    c3(aa)=16/(9*(1+nx(aa)));
    m_alpha(aa,1)=f*(u0*Ms)^2*(6*nx(aa)*average4_R/(8*average3_R)-2/3*c2*mx(aa)^2*average3_R^(1/3));
    m_alpha(aa,2)=f*(u0*Ms)^2*(3*((3/4*nny(aa)+1/4*nnz(aa))/numb)*average4_R/(8*average3_R)-2/3*c2*my(aa)^2*average3_R^(1/3));
    m_alpha(aa,3)=f*(u0*Ms)^2*(3*((1/4*nny(aa)+3/4*nnz(aa))/numb)*average4_R/(8*average3_R)-2/3*c2*mz(aa)^2*average3_R^(1/3));
    
    %kesa(aa)=2*f*(u0*Ms)^2*average4_R/(3*average3_R*c3(aa))*(1-c2*c3(aa)*mm(aa)^2*average3_R^(4/3)/average4_R);
    kesa0(aa)=m_alpha(aa,1)+m_alpha(aa,2)+m_alpha(aa,3);
    md_cos0(aa,1)=m_alpha(aa,1)/kesa0(aa);% mean x-direction cosione, coincide with gamma in Suzanne's work
    md_cos0(aa,2)=m_alpha(aa,2)/kesa0(aa);% mean y-direction cosione, coincide with gamma in Suzanne's work
    md_cos0(aa,3)=m_alpha(aa,3)/kesa0(aa);% mean z-direction cosione, coincide with gamma in Suzanne's work
    nd_angle0(aa)=f*mm(aa)*Bs*L*sqrt(c)/1e6; % Rotation angle of polarization [rad]
    nd_delta0(aa)=3*c3(aa)*kesa0(aa)/(2*f*(u0*Ms)^2*(1-c2*c3(aa)*mm(aa)^2));%[um]
end

%{
if aax0>1
  if direction==1
     figure
     plot(Ba(aa,3),mz(aa),'ro');
     xlabel('Applied field Ba (T)','FontSize',20);
     ylabel('Reduced magnetization','FontSize',20);
     hold on;
  else
     plot(Ba(aa,3),mz(aa),'b+');
     xlabel('Applied field Ba (T)','FontSize',20);
     ylabel('Reduced magnetization','FontSize',20);
     title('Magnetic hysteresis','FontSize',20);
     hold on;
  end
end
%}


Mag_parameter0(aa,:)=[Ba(aa,3);Ba(aa,3)/(u0*Ms);mx(aa);my(aa);mz(aa);mm(aa);nx(aa);ny(aa);nz(aa);M_total(cc-1,1);M_total(cc-1,2);M_total(cc-1,3)];
%xlswrite('D:\hfang\Desktop\results.xlsx',Mag_parameter,'Magnetic','B1:B12');% save the matrix of magnetization
%xlswrite('D:\hfang\Desktop\results.xlsx',Mk,'Magnetic','C2:E31');% save the matrix of magnetization
%}
aa
end

% Save outputs
% if direction==1
%    xlswrite('D:\hfang\Desktop\results.xlsx',Mag_parameter,'Magnetic','H2:S41');% save the matrix of magnetization
%    xlswrite('D:\hfang\Desktop\results.xlsx',Mk,'Magnetic','u2:w301');% save the matrix of magnetization
% else
%    xlswrite('D:\hfang\Desktop\results.xlsx',Mag_parameter,'Magnetic','AA2:AL41');% save the matrix of magnetization
%    xlswrite('D:\hfang\Desktop\results.xlsx',Mk,'Magnetic','AN2:AP301');% save the matrix of magnetization
% end

if aax0>1
   if aa==aax0
      direction=-direction;
      stop_flag=stop_flag-1;
   end
else
    if aa==aax0
       stop_flag=stop_flag-2;
    end
end
% xlswrite('D:\hfang\Desktop\results.xlsx',E_total,'Sheet5','A1:Ad21');
% xlswrite('D:\hfang\Desktop\results.xlsx',fai_solution,'Sheet6','A1:Ad21');
end
end

%{
Mag_anist=[4 2 1];
Hk=[15 97 42];
Ms=1e4;
fai=pi/3;
    f=@(N)magnetization_vector(N,Mag_anist,Hk,Ms,fai);
    [out,fval,exitflag]=fsolve(f,[3 7]);
    Mk(1)=out(1)*Mag_anist(1)+out(2)*Hk(1);% Updated magnetization for each particle in X direction
    Mk(2)=out(1)*Mag_anist(2)+out(2)*Hk(2);% Updated magnetization for each particle in Y direction
    Mk(3)=out(1)*Mag_anist(3)+out(2)*Hk(3);% Updated magnetization for each particle in Z direction
sqrt(Mk(:,1).^2+Mk(:,2).^2+Mk(:,3).^2)
out
Mk
%}

%{
for jj=1:numb;
    plot3(Mag_anist(jj,1),Mag_anist(jj,2),Mag_anist(jj,3),'ro');
    title('Vector of the orientation');
    hold on;
end
%}
