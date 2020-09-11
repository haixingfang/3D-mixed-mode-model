% This is the main function: ferrite_3d_model_voronoi_PBC_ND_CNT.m
% voronoi structure, PBC,spherical ferrite, Classical nucleation theory,
% mixed-mode, para-equilibrium, effective interface mobility
% coupling to magnetic and ND models
% Copyright by Haixing Fang, Nov 2016
% At present: only valid for Fe-xC-0.49Mn (wt.%)
% Use Fe-0.1C-0.49Mn (wt.%) as an example

clear all;
close all;
% load Multi-Parametric Toolbox 3
% https://www.mpt3.org/
addpath('C:\Users\hfang\Documents\MATLAB\tbxmanager');
if isdir('.\tbxmanager')
    addpath('.\tbxmanager');
    tbxmanager restorepath;
    mpt_init;
else
    fprintf(['The path of mpt3 toolbox is not found. \n', ...
        'Please install mpt3 toolbox first! \n', ...
        'ttps://www.mpt3.org/ \n']);
end
% cite using MPT3:
% M. Herceg, M. Kvasnica, C.N. Jones, and M. Morari. Multi-Parametric Toolbox 3.0. In Proc. of the European Control Conference, pages 502?10, Zurich, Switzerland, July 17?9 2013.

% characteristics of Fe-0.1C-0.49Mn
T_a3=1106;   % A3-temperature othorequilibrium=1125K, Para-equilibrium=1117K, linerization of PE from M.G.Meccozi=1106K [K]
T_a1=987; % A1-temperature othorequilibrium [K]
Tc=1043-10*0.49;     % Curie temperature [K]
R=8.314; % gas constant [J/(mol.K)]

TR=1073; % reference temperature [K] 
wCR_bcc=0.009; % wC in bcc at TR [wt.%]
wCR_fcc=0.279; % wC in fcc at TR [wt.%]
m_bcc=-10350;  % slope of linear fits for (fcc+bcc)/bcc boundary [K/wt.%]
m_fcc=-186.2;  % slope of linear fits for fcc/(fcc+bcc) boundary[K/wt.%]
deltaS=3.5e5;  % Entropy [J/(K.m3)]
T_wC=[75.12902 -268.19506 1143.7567]; % TA3=a*wC^2+b*wC+c [K]
% T_a3=T_wC(1)*Comp(1)^2+T_wC(2)*Comp(1)+T_wC(3); % A3 [K]

% modelling parameters
T_a3=1116; % [K]
T0=T_a3; % starting temperature [K]
rou_bcc=0.6e14; % ferrite nuclei density [m-3]
% Rate=10;  % cooling rate [K/s]
Rate=0.4;  % cooling rate [K/s]
u=0.13;   % to adjust the growth velocity
eps=1e-6; % minimum error

% starting microstrucuture
Lb=70; % box size [um]
Npot=round(rou_bcc*(Lb*1e-6)^3); % maximum number of ferrite nuclei
f_N=1;       % scaling factor
N0=Npot*f_N; % scaled number of potential nucleation sites [-]
dmin=12;     % The minimum distance required between the centroids of austenite grains
dfcc=20; % average austenite grain diameter [um]
N_fcc=round(Lb^3/(pi*(dfcc^3)/6)); % number of austenite grains
rou_fcc=N_fcc/Lb^3; % Number density of austenite grain [um-3]
shieldD=1/4*(Lb^3/N0)^(1/3); % 1/4 of the everage ferrite spacing [um]

enlarge_Lb=Lb+2*(0.5*1/rou_fcc^(1/3));% enlarged by 2*average neibouring distance
enlarge_N=round(rou_fcc*enlarge_Lb^3); % Number of austenite grains in enlargement box (enlarged by 2*average neibouring distance) that have the same number density
minus_edge=0.5*1/rou_fcc^(1/3);

% Define the anistropy axis for ferrite nucleation sites
MA_theta0  = 2*pi*rand(1,Npot);
MA_phi0    = pi*rand(1,Npot);

% Assign the initial value for magnetization of each ferrite grain
alpha0=2*pi*rand(Npot,1);
beta0=pi*rand(Npot,1);

%General constantants
M_Fe=56;
M_C=12;
M_Mn=55;
M_Si=28;
%Comp=[0.364 0.656 0.305];%C Mn Si in wt.% the same as in Suzanne'work for C35 Steel
Comp=[0.1 0.49 0];%C Mn Si in wt.% the same as in M.Militzer Acta Mater. 2006 Phase field modelling of austenite to ferrite transformation
% Comp=[0.247 2.06 0.098];%C Mn Si in wt.% of the sample collected from Ancelor Mittar
Mole=Comp(1)/M_C+Comp(2)/M_Mn+Comp(3)/M_Si+(100-Comp(1)-Comp(2)-Comp(3))/M_Fe;%total mole number of elements
Comp_m(1)=Comp(1)/M_C/Mole;Comp_m(2)=Comp(2)/M_Mn/Mole;Comp_m(3)=Comp(3)/M_Si/Mole;%mole fraction of each element
Comp_m=100.*Comp_m; % [mol%]
Ux=Comp_m(2)/(Comp_m(2)+(100-Comp_m(1)-Comp_m(2)-Comp_m(3))); % U-fraction X/(Fe+X)

%coefficents for fitting equilibrium lines calculated from TC
%line1:Fe-C;line2:Fe-C-0.656Mn;line3:Fe-C-0.656Mn-0.305Si(Suzanne);line4:Fe-C-0.611Mn(Offerman);line5:
%Fe-C-0.49Mn (M.Militzer);line6: Fe-C-2.06Mn(Ancelor Mittar)
a_eq=[1.42754e-5 -0.03519 21.66617;0.00001871 -0.04471291 26.70387911;0.000014365 -0.03537227 21.73769738;0.00001632 -0.03956033 23.94207419; ...
    0.000015666 -0.038180593 23.221586256;0.000027216 -0.062554358 35.892301031];%A3-line,binomial fitting,f(x)=ax^2+bx+c,
b_eq=[-0.000118 0.139638;-0.00006564 0.07518773;-0.00007891 0.09090873;-0.00004104 0.04897501; ...
    -0.00008864 0.10182309;-0.00003886 0.03977055];%carbon solubility in ferrite,linear fitting,f(x)=ax+b,
Temp_eq=[1068 1000;1053 984;1064 985;1054 985;1125 986;1045 966];%colume1:A3-temperature;colume2:A1+temperature,line has the same meaning as above matrix

%%%%%%%%%%%%% phase diagram in para-equilibrium condition Nov21, 2016
a_eq(5,:)=[0.0000166618 -0.04012 24.12059];
b_eq(5,:)=[-0.0000996012 0.11565];

%%% 2-nd polynomial fits for kafang as a function of T [K]
kafang_p=[0.0025 -5.365 2979.3]; % kafang(T)=a*T^2+b*T+c [J/(mol*mol%)]
% Meff=3.5e-7; % effective interface mobility [m3.m/(J.s)] in Pina and C.Bos Meff=2.5~5e-7, E.Gamsjager think Meff=6.9~10e-7
M0=2.4e-7; % [m3.m/(J.s)]
Meff=M0/7.1e-6; % [mol.m/(J.s)]
QM=140e3; % [J/mol]

%Generate positions potential nucleation sites and plot the 3d voronoi diagram
generate_new_seed=0; % generate new seeding points: 1-yes; 0-no; by defaul to load the old seeding points
if generate_new_seed==1
    A1=enlarge_Lb*rand([1,3])-minus_edge;
    i=1;
    A=A1;
    while i<=enlarge_N-1
       A2=enlarge_Lb*rand([1,3])-minus_edge-Lb/2;
       comb=[A;A2];
       distance=pdist(comb,'euclidean');
       if min(distance)>dmin
          A=[A;A2];
          i=i+1
       end
    end
    % export the predefined austenite centroid coordinates
    dlmwrite(strcat(num2str(length(A)),'AusteniteCentroidsSymmetric.txt'),[A(:,1) A(:,2) A(:,3)],'delimiter',' ');
else
    % load the predefined austenite centroid coordinates
    FileName_prefix='153AusteniteCentroidsSymmetric'; % randomly generate 153 symmetric centers
    FileName_pattern='.txt';
    baseFileName=[FileName_prefix FileName_pattern];
    fullFileName=fullfile(pwd, baseFileName);
    fileID=fopen(fullFileName,'r');
    A=[];
    while(~feof(fileID))
        textdata=str2num(fgetl(fileID));
        if isempty(textdata)
            continue;
        else
            A=[A; textdata];
        end
    end
    fclose(fileID);
end
dis=pdist(A,'euclidean'); % Pairwise distance of austenite grain center
B=Polyhedron([-minus_edge-Lb/2 -minus_edge-Lb/2 -minus_edge-Lb/2; ...
    enlarge_Lb-minus_edge-Lb/2 -minus_edge-Lb/2 -minus_edge-Lb/2; ...
    enlarge_Lb-minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2 -minus_edge-Lb/2; ...
    -minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2 -minus_edge-Lb/2; ...
    -minus_edge-Lb/2 -minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2; ...
    enlarge_Lb-minus_edge-Lb/2 -minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2; ...
    enlarge_Lb-minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2; ...
    -minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2]);% Boundary vertices
[V,cells] = mpt_voronoi(A', 'bound', B); % Call for the mpt_voronoi
vertices = V.Set.forEach(@(e) e.V, 'UniformOutput', false);% find the vertices for each cell
new_vertices = cat(1, vertices{:}); % Combine the vertices into one matrix
potential0=unique(round(new_vertices*1e6), 'rows')/1e6; % Get rid of the numerical noise 
to_remove = [-Lb/2 -Lb/2 -Lb/2;Lb/2 -Lb/2 -Lb/2;Lb/2 Lb/2 -Lb/2;-Lb/2 Lb/2 -Lb/2; ...
    -Lb/2 -Lb/2 Lb/2;Lb/2 -Lb/2 Lb/2;Lb/2 Lb/2 Lb/2;-Lb/2 Lb/2 Lb/2];% Remove the corners of the Lb*Lb*Lb cubic
potential=setdiff(potential0, to_remove, 'rows');
for i=1:length(potential)
    for j=1:3
    if potential(i,j)>Lb/2
       potential(i,j)=Lb/2+1;
    else if potential(i,j)<-Lb/2
        potential(i,j)=-Lb/2-1;
        end
    end
    end
end
site=[];
% Remove points beyond the Lb*Lb*Lb cubic box
for i=1:length(potential)
  if all(potential(i,:)-Lb/2-1)&&all(potential(i,:)+Lb/2+1)
     eval(['site',num2str(i),'=','potential(i,:)']); % Convert the number to string
     eval(['site=[site;site',num2str(i),'];']);  % Combine site1, site 2,...
  end
end
p1=randperm(length(site));
position=site(p1(1:N0),:);% Randomly select N0 potential nucleation site for ferrite
site=site(randperm(length(site)),:); % randomly sort the possible nucleation sites

figure('Name','Austenite geometry');
subplot(1,2,1);
V.plot('alpha', 0.13); % Adjust the transparancy
hold all;
plot3(A(:,1),A(:,2),A(:,3),'b+','LineWidth',2);% Plot the centroids
plot3(potential0(:,1),potential0(:,2),potential0(:,3),'r.','MarkerSize',20); % Plot all the corners
axis([-enlarge_Lb/2 enlarge_Lb/2 -enlarge_Lb/2 enlarge_Lb/2 -enlarge_Lb/2 enlarge_Lb/2]);
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
zlabel('z (\mum)','FontSize',16);
set(get(gca,'xlabel'),'rotation',18);
set(get(gca,'ylabel'),'rotation',-25);
set(get(gca,'zlabel'),'rotation',90);
set(gca,'fontsize',14);
set(gca,'linewidth',2);
% legend('+ centroids of austenite grains','. corners of austenite grains');
box on;
grid off;
hold off;

subplot(1,2,2);
V.plot('alpha', 0.13); % Adjust the transparancy
hold all;
plot3(A(:,1),A(:,2),A(:,3),'b+','LineWidth',2);% Plot the centroids
plot3(site(:,1),site(:,2),site(:,3),'r.','MarkerSize',20);% Plot the corners
% plot3(position(:,1),position(:,2),position(:,3),'b*','LineWidth',2);
grid off;
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
zlabel('z (\mum)','FontSize',16);
axis([-Lb/2 Lb/2 -Lb/2 Lb/2 -Lb/2 Lb/2]);
% axesLabelsAlign3D;
set(get(gca,'xlabel'),'rotation',18);
set(get(gca,'ylabel'),'rotation',-25);
set(get(gca,'zlabel'),'rotation',90);
set(gca,'fontsize',14);
set(gca,'linewidth',2);
box on;
legend('+ centroids of austenite grains','. potential ferrite nucleation sites','* randomly selected nucleation sites');
hold off;

% Define initial features for each austenite
for i=1:length(A)
   A_P(i,1)=i; % order number
   A_P(i,2)=V.Set(i).Data.voronoi.seed(1,1); % x coordinate
   A_P(i,3)=V.Set(i).Data.voronoi.seed(2,1); % y coordinate
   A_P(i,4)=V.Set(i).Data.voronoi.seed(3,1); % z coordinate
   A_P(i,5)=V.Set(i).volume; % volume [um3]
   A_P(i,6)=0; % number of nucleated ferrite in its vertex
   A_P(i,7)=length(V.Set(i).V); % total number of vertex
   A_P(i,8)=0;
   for j=1:length(V.Set(i).V)
       if abs(V.Set(i).V(j,1))>Lb/2 || abs(V.Set(i).V(j,2))>Lb/2 || abs(V.Set(i).V(j,3))>Lb/2
          A_P(i,8)=A_P(i,8)+1; % total number of vertex outside Lb
       end
   end
   A_P(i,9)=A_P(i,7)-A_P(i,8); % effective number of vertex(<Lb)
   A_P(i,10)=Comp(1); % C content in matrix [wt.%]
   A_P(i,11)=Comp(2); % Mn content in matrix [wt.%]
   A_P(i,12)=0; % ferrite volume in this specific austenite [um3]
   A_P(i,13)=A_P(i,12)/A_P(i,5); % ferrite fraction of each austenite
   A_P(i,14)=0; % soft impingement_1 or not_0
   A_P(i,15)=Comp_m(1); % C content in matrix [mol%]
   A_P(i,16)=Comp_m(2); % C content in matrix [mol%]
end
[position_MXY position_MXZ position_MYZ position_MO]=MirrorPoints(site,Lb); % Get the coordinates of the mirror points

%Generate a N0*10 matrix:no,x,y,z,active or not,nucleation
%time,radius,effective radius,nucleated flag,actual volume,ratio of actual
%volume to the original volume
for i=1:length(site)
    N_p(i,1)=i;%site number
    N_p(i,2)=site(i,1);%x coordinate
    N_p(i,3)=site(i,2);%y coordinate
    N_p(i,4)=site(i,3);%z coordinate
    N_p(i,5)=1; % active=1, not active=0
    N_p(i,6)=-1;% nucleation time, not nucleated=-1
    N_p(i,7)=0; % original radius
    N_p(i,8)=0; % effective radius
    N_p(i,9)=0; % not nucleated=0, nucleated=1
    N_p(i,10)=0;% actual volume (after substract the overlay volume)
    N_p(i,11)=-1;% actual volume/original volume,=-1 when original volume=0
    N_p(i,12)=0;% actual volume fraction/equilibrium volume fraction predicted by phase diagram
    N_p(i,13)=0;% The flag of impingement of 4 spheres (0 denotes no impingement while 1 denotes impingement)
    N_PR{i}=[];
    for j=1:length(A)
      for k=1:length(V.Set(j).V)
       if abs(site(i,1)-V.Set(j).V(k,1))<=eps && abs(site(i,2)-V.Set(j).V(k,2))<=eps && abs(site(i,3)-V.Set(j).V(k,3))<=eps
          N_PR{i}=[N_PR{i}; j k]; % derive the corresponding relationship between vertex and austenite grains
       end
      end
    end
    for m=1:length(N_PR{i}(:,1))
        N_PR{i}(m,3)=sqrt((site(i,1)-A_P(N_PR{i}(m,1),2))^2+(site(i,2)-A_P(N_PR{i}(m,1),3))^2+(site(i,3)-A_P(N_PR{i}(m,1),4))^2);% identical distance to neighbouring centers
    end
    N_p(i,14)=mean(N_PR{i}(:,3)); % average distance between vertex and austenite centers [um]
    N_p(i,15)=mean(A_P(N_PR{i}(:,1),10)); % average remote C content in austenite [wt.%]
    N_p(i,16)=mean(A_P(N_PR{i}(:,1),11)); % average remote Mn content in austenite [wt.%]
    N_p(i,17)=0.004; % average remote C content in ferrite at T_a3 [wt.%]
    N_p(i,18)=mean(A_P(N_PR{i}(:,1),11)); % average remote Mn content in ferrite [wt.%]
    N_p(i,19)=T0-T_a3; % local undercooling [K]
    N_p(i,20)=0; % G_fcc-G_bcc
    N_p(i,21)=1; % only becomes 1 means it can be potential nucleation sites
    N_p(i,22)=Comp(1); % C content at the site [wt.%]
    N_p(i,23)=Comp(2); % Mn content at the site [wt.%]
    N_p(i,24)=0; % interface velocity [um/s]
    N_p(i,25)=0; % flag of soft impingement either 0 or 1
    N_p(i,26)=0; % hard impingement either 0 or 1
    N_p(i,27)=N_p(i,14); % initial soft impingement distance [um]
    if i==1
       N_p(i,28)=min([mean(pdist([site(i,:);site(i+1:end,:)])) mean(pdist([position_MXY(i,:);site(i+1:end,:)])) ...
           mean(pdist([position_MXZ(i,:);site(i+1:end,:)])) mean(pdist([position_MYZ(i,:);site(i+1:end,:)])) ...
           mean(pdist([position_MO(i,:);site(i+1:end,:)]))]); % mean distance between vertices [um]
    else if i==length(site)
       N_p(i,28)=min([mean(pdist([site(i,:);site(1:i-1,:)])) mean(pdist([position_MXY(i,:);site(1:i-1,:)])) ...
           mean(pdist([position_MXZ(i,:);site(1:i-1,:)])) mean(pdist([position_MYZ(i,:);site(1:i-1,:)])) ...
           mean(pdist([position_MO(i,:);site(1:i-1,:)]))]); % mean distance between vertices [um]
        else
       N_p(i,28)=min([mean(pdist([site(i,:);site(1:i-1,:)])) mean(pdist([position_MXY(i,:);site(1:i-1,:)])) ...
           mean(pdist([position_MXZ(i,:);site(1:i-1,:)])) mean(pdist([position_MYZ(i,:);site(1:i-1,:)])) ...
           mean(pdist([position_MO(i,:);site(1:i-1,:)]))]); % mean distance between vertices [um]
       N_p(i,28)=mean([N_p(i,28) min(pdist([site(i,:);site(i+1:end,:)])) min(pdist([position_MXY(i,:);site(i+1:end,:)])) ...
           min(pdist([position_MXZ(i,:);site(i+1:end,:)])) min(pdist([position_MYZ(i,:);site(i+1:end,:)])) ...
           min(pdist([position_MO(i,:);site(i+1:end,:)]))]); % minimum distance between vertices [um]
        end
    end
    N_p(i,29)=0; 
    for m=1:length(N_PR{i}(:,1))
        N_p(i,29)=N_p(i,29)+A_P(N_PR{i}(m,1),5); % volume that it can maximum grow [um^3]
    end
    N_p(i,30)=-1; % restore the moment when extended volume is corrected [s]
    N_p(i,31)=0;  % restore the site number of pair growing ferrite
    N_p(i,32)=(3.*sum(A_P(N_PR{i}(:,1),5))/(4*pi)).^(1/3); % the radius of remaining untransformed surrounding austenite volume [um]
    N_p(i,33)=0;
    N_pD(i,7)=0; % define an initial value for grain radius in diffusion-controlled mode [um]
    i
end
for j=1:length(A)
    A_PR{j}=[];
    for k=1:length(V.Set(j).V)
        for i=1:length(N_p(:,1))
           if abs(N_p(i,2)-V.Set(j).V(k,1))<=eps && abs(N_p(i,3)-V.Set(j).V(k,2))<=eps && abs(N_p(i,4)-V.Set(j).V(k,3))<=eps
              A_PR{j}(k)=N_p(i,1); % restore the site order number
           end
        end
    end
end

Timer(1)=0; % total crystallization time
i=1;j=1;k=1; % initial definition
t=0;    % define initial t
dT=1; dt=1*dT/Rate; % T and time step of calculation in seconds
start=0; % initiation of the first nuclei time
dd=1; % Count number for the calculation of magnetic_m
stop_cycle=0;
F(1)=0; % initial ferrite fraction
sumN(1)=0; % total ferrite nucleis
MagFlag=0; % whether to calculate ND characters: 0-no; 1-yes. By default it is set as 0.

% main loop to calculate nucleation,growth and call impingement function
while stop_cycle~=1 
    t
    n=1;
    Nucleated{i}=[]; % restore all the nucleated ferrite
    N_pA{i}=[];
    Timer(i)=t; % timer [s]
    T(i)=T0-Rate*t;%temperature at time t
    for m=1:length(N_p(:,1))
%         N_p(m,19)=T(i)-(TR+0.5*(m_fcc*(N_p(m,15)-wCR_fcc)+m_bcc*(N_p(m,17)-wCR_bcc))); % local undercooling [K], Pina's method
        N_p(m,19)=T(i)-(T_wC(1)*N_p(m,22)^2+T_wC(2)*N_p(m,22)+T_wC(3));
        N_p(m,20)=-deltaS*N_p(m,19);% driving force [J/m3] positive means enabling austenite to ferrite
        if N_p(m,20)<0
            N_p(m,21)=0; % only becomes 1 means it can be potential nucleation sites
        else
            N_p(m,21)=1;
        end
        if N_p(m,9)==0 && N_p(m,5)==1
           N_pA{i}(n,:)=N_p(m,:); % restore the potential available site
           n=n+1;
        end
    end
    if ~isempty(N_pA{i})
        deltaGV(i)=mean(N_pA{i}(:,20)); % driving force averaging the whole structure [J/m3]
    else
        deltaGV(i)=0;
    end
    if Timer(i)==0
        Nt(i)=0;
    else
%         Nt(i)=CNT(T(i),F(i-1),sum(N_p(:,5))-sum(N_p(:,9)),deltaGV(i),Timer(i)); % nucleation rate based on Classical Nucleation Theory [s-1]
        Nt(i)=PhaseFieldNucleation(T(i),Rate,dt,N0); % nucleation same as M.G.Meccozi used in phase field approach [s-1]
    end
    N(i)=Nt(i)*dt; % nucleation number in dt
    N(i)=round(N(i))
    if i>=2
        sumN(i)=sumN(i-1)+N(i);
        if sumN(i)>=Npot
            sumN(i)=Npot;
            N(i)=sumN(i)-sumN(i-1);
        end
    end

    N_eff(i)=0;%N_eff is the real nucleation number
    if i>1 && sumN(i-1)==0 && sumN(i)>0
        start=Timer(i);%start moment of the 1st nuclei
    end
    numb=1;
    if ~isempty(N_pA{i})
    [sortN_p{i},sortIndex{i}] = sortrows(N_pA{i}, 22);% min wC should nucleat first
    if t>=start
        j=1;
       while (isempty(Nucleated{i})&&N(i)>0) || (~isempty(Nucleated{i})&&length(Nucleated{i}(:,1))<=N(i)-1)
           Index=N_pA{i}(sortIndex{i}(j),1); % site oder number
           % Update the distance from the nuclei
           if ~isempty(Nucleated{i})
               N_p(Index,28)=min([N_p(Index,28) min(pdist([N_p(Index,2:4);Nucleated{i}(1:end,2:4)]))]);
           end
           if N_p(Index,5)==1 && N_p(Index,21)==1 && N_p(Index,9)==0 && N_p(Index,28)>shieldD% active state and driving force > 0 and not nucleated
             N_p(Index,6)=Timer(i);%the specific nucleation time for particle j
             N_p(Index,7)=5e-3; % nucleus radius [um]
             N_p(Index,8)=5e-3; % nucleus radius [um]
             N_pD(Index,7)=5e-3; % nucleus radius [um] for diffusion-controlled mode
             N_pD(Index,8)=5e-3; % nucleus radius [um] for diffusion-controlled mode
             N_p(Index,9)=1;% particle j nucleated
             Nucleated{i}=[Nucleated{i};N_p(Index,:)]; % restore the nucleated ferrite
             N_eff(i)=N_eff(i)+1;%N_eff is the real nucleation number
             j=j+1;
           else
              j=j+1;
           end
           if j>length(sortIndex{i})
               break;
           end
       end
       if i>1
          N_eff(i)=N_eff(i)+N_eff(i-1);
          Nucleated{i}=[Nucleated{i-1};Nucleated{i}]; % combine all the nucleated ferrite
       end
    end
    else
        Nucleated{i}=[Nucleated{i-1};Nucleated{i}];
        N_eff(i)=N_eff(i-1);
    end
    
    if ~isempty(Nucleated{i})
       numb=length(Nucleated{i}(:,1));% let numb=j,reveals the number of nucleated particle
       if numb>N0
          numb=N0;
       end
    else
        numb=1;
    end

    if T(i)<T_a3  %When T is above A3-temperature, there is no C redistribution
        wC_A_eq(i)=a_eq(5,1)*T(i)^2+a_eq(5,2)*T(i)+a_eq(5,3);% This equation of A3-line is fitted from TC data
        wC_F_eq(i)=b_eq(5,1)*T(i)+b_eq(5,2);% This equation of C solubility in ferrite is fitted from TC data
    else
        wC_A_eq(i)=Comp(1);
        wC_F_eq(i)=0;
    end
    xC_A_eq(i)=100*(wC_A_eq(i)/M_C)/(wC_A_eq(i)/M_C+(100-wC_A_eq(i))*Ux/M_Mn+(100-wC_A_eq(i))*(1-Ux)/M_Fe); % [mol%]
    xC_F_eq(i)=100*(wC_F_eq(i)/M_C)/(wC_F_eq(i)/M_C+(100-wC_F_eq(i))*Ux/M_Mn+(100-wC_F_eq(i))*(1-Ux)/M_Fe); % [mol%]
    F_eq(i)=(wC_A_eq(i)-Comp(1))/(wC_A_eq(i)-wC_F_eq(i));%Equilibrium ferrite fraction predicated by the phase diagram at time t
    if F_eq(i)<0
        F_eq(i)=1e-12;
    end

    % below is to calculate the ferrite growth
    if i==1
       x_C=Comp_m(1)/100; % x_C in at.
    else
       x_C=(Comp_m(1)-F(i-1)*xC_F_eq(i-1))/(1-F(i-1))/100; % x_C in at.
    end
    y_C=x_C/(1-x_C);
    
    D_C(i)=4.53e-7*((1+y_C*(1-y_C)*8339.9/T(i))*exp(-(1/T(i)-2.221e-4)*(17767-26436*y_C)))*1e12;% Volume diffusion of Carbon in um2/s, J Agren 1986
    D_C1(i)=2.343e-5*exp(-148e3/(R*T(i)))*1e12; % Volume diffusion of Carbon in austenite um2/s, J Agren 1986
    D_C2(i)=1.5e-5*exp(-142.1e3/(R*T(i)))*1e12; % Volume diffusion of Carbon in austenite um2/s, R.C.Weast, 1989; C.Bos and J. Sietsma, 2007
    Mobility(i)=Meff*exp(-QM/(R*T(i)));% effective interface mobility [mol.m/(J.s)]
    Kafang(i)=kafang_p(1)*T(i)^2+kafang_p(2)*T(i)+kafang_p(3); % kafang(T)=a*T^2+b*T+c [J/(mol*mol%)]
%     plot(T,D_C,'ro-',T,D_C1,'b.-',T,D_C2,'m*-')
%     legend('site fraction-J Agren 1986','Arrhenius-J Agren 1986','Arrhenius-R C Weast 1989')

    ll=1;
    if i>1 && ~isempty(Nucleated{i})
    while ll<=length(Nucleated{i}(:,1))
      l=Nucleated{i}(ll,1); % site oder number
      if F(i-1)<=1.01*F_eq(i)
         if t>=N_p(l,6) && N_p(l,6)>0 % && N_p(l,13)==0
             for r=1:length(N_PR{l}(:,1))
                 if N_p(l,26)==0  % no hard impingement
                     [DiffInfo{l}(i,r,:)]=Mixmode_diffusion_profile(xC_F_eq(i),xC_A_eq(i),Comp_m(1),D_C(i)*1e-12,Mobility(i),Kafang(i),N_p(l,7),length(N_PR{l}(:,1)),N_p(l,32));% mix-mode
                     HardFlag{l}(i,r)=0; % hard flag
                     [DiffInfoD{l}(i,r,:)]=Diffusion_control_profile(xC_F_eq(i),xC_A_eq(i),Comp_m(1),D_C(i)*1e-12,Mobility(i),Kafang(i),N_pD(l,7),length(N_PR{l}(:,1)),N_p(l,32)); % diffusion-controlled
                     HardFlagD{l}(i,r)=0; % hard flag
                     if DiffInfo{l}(i,r,2)<0 || DiffInfo{l}(i,r,2)>xC_A_eq(i) % fail to identify the hard impingement
                         DiffInfo{l}(i,r,1)=(Comp_m(1)-N_p(l,12)*xC_F_eq(i))/(1-N_p(l,12));
                         DiffInfo{l}(i,r,2)=DiffInfo{l}(i,r,1);
                         DiffInfo{l}(i,r,3)=0;
                         DiffInfo{l}(i,r,7)=1;
                         N_p(l,25)=1;
                         N_p(l,26)=1;
                     end
                     if DiffInfo{l}(i,r,2)<DiffInfo{l}(i-1,r,2)
                         DiffInfo{l}(i,r,2)=DiffInfo{l}(i-1,r,2)+eps/2;
                     end
                 
                     MeanField=(Comp_m(1)-N_p(l,12)*xC_F_eq(i))/(1-N_p(l,12));
                     if DiffInfo{l}(i,r,2)<=DiffInfo{l}(i,r,1)
                         DiffInfo{l}(i,r,2)=MeanField;
                         DiffInfo{l}(i,r,1)=MeanField;
                         DiffInfo{l}(i,r,3)=0;
                         DiffInfo{l}(i,r,4)=N_p(l,7);
                         DiffInfo{l}(i,r,5)=xC_F_eq(i);
                         DiffInfo{l}(i,r,6)=xC_A_eq(i);
                         DiffInfo{l}(i,r,7)=1;
                         HardFlag{l}(i,r)=1; % hard flag
                         N_p(l,33)=1; % mean field flag
                     end
                     if DiffInfo{l}(i,r,3)<0
                         DiffInfo{l}(i,r,3)=N_p(l,27)-DiffInfo{l}(i,r,4);
                     end
                     if DiffInfoD{l}(i,r,2)<0 || DiffInfo{l}(i,r,2)>xC_A_eq(i)
                         DiffInfoD{l}(i,r,2)=xC_A_eq(i);
                     end
                     if DiffInfoD{l}(i,r,3)<0
                         DiffInfoD{l}(i,r,3)=N_pD(l,27)-DiffInfoD{l}(i,r,4);
                     end
                     % DiffInfo=[Xpm Xip DiffLL Rbcc Xneq Xpeq softflag];
                     % mol%, mol%, um, um, mol%, mol%, 1or0;
                 else
                     DiffInfo{l}(i,r,1)=(Comp_m(1)-N_p(l,12)*xC_F_eq(i))/(1-N_p(l,12));
                     if DiffInfo{l}(i,r,1)>xC_A_eq(i)
                         DiffInfo{l}(i,r,1)=xC_A_eq(i); % stop when reach equilibrium
                     else if DiffInfo{l}(i,r,1)<DiffInfo{l}(i-1,r,1)
                             DiffInfo{l}(i,r,1)=DiffInfo{l}(i-1,r,1)+eps/2;
                         end
                     end
                     DiffInfo{l}(i,r,2)=DiffInfo{l}(i,r,1); % Xip [mol%] 
                    if DiffInfo{l}(i,r,2)<DiffInfo{l}(i-1,r,2)
                        [DiffInfo{l}(i,r,:)]=Mixmode_diffusion_profile(xC_F_eq(i),xC_A_eq(i),A_P(N_PR{l}(r,1),15),D_C(i)*1e-12,Mobility(i),Kafang(i),N_p(l,7),length(N_PR{l}(:,1)),N_p(l,32));% mix-mode
                        if DiffInfo{l}(i,r,1)>xC_A_eq(i)
                            DiffInfo{l}(i,r,1)=xC_A_eq(i);
                        end
                        if DiffInfo{l}(i,r,2)>xC_A_eq(i)
                            DiffInfo{l}(i,r,2)=DiffInfo{l}(i,r,1);
                        end
                    else
                        DiffInfo{l}(i,r,3)=0; % DiffL [um]
                        DiffInfo{l}(i,r,4)=N_p(l,7); % original radius [um]
                        DiffInfo{l}(i,r,5)=xC_F_eq(i); % Xneq [mol%]
                        DiffInfo{l}(i,r,6)=xC_A_eq(i); % Xpeq [mol%] 
                    end   
                     DiffInfo{l}(i,r,7)=1; % soft flag
                     HardFlag{l}(i,r)=1; % hard flag
                 
                     DiffInfoD{l}(i,r,1)=(Comp_m(1)-A_P(N_PR{l}(r,1),13)*xC_F_eq(i))/(1-A_P(N_PR{l}(r,1),13)); % Xpm [mol%]
                     if DiffInfoD{l}(i,r,1)>xC_A_eq(i)
                         DiffInfoD{l}(i,r,1)=xC_A_eq(i); % stop when reach equilibrium
                         else if DiffInfoD{l}(i,r,1)<DiffInfoD{l}(i-1,r,1)
                             DiffInfoD{l}(i,r,1)=DiffInfoD{l}(i-1,r,1);
                             end
                     end
                     DiffInfoD{l}(i,r,2)=xC_A_eq(i); % Xip [mol%]
                     DiffInfoD{l}(i,r,3)=0; % DiffL [um]
                     DiffInfoD{l}(i,r,4)=N_pD(l,7); % original radius [um]
                     DiffInfoD{l}(i,r,5)=xC_F_eq(i); % Xneq [mol%]
                     DiffInfoD{l}(i,r,6)=xC_A_eq(i); % Xpeq [mol%]
                     DiffInfoD{l}(i,r,7)=1; % soft flag
                     HardFlagD{l}(i,r)=1; % hard flag
                 end
             end
            v_t(i,l)=Mobility(i)*Kafang(i)*(xC_A_eq(i)-mean(DiffInfo{l}(i,:,2)))*1e6; % velocity [um/s] mix-mode
            if v_t(i,l)<0
                v_t(i,l)=0;
            end
            N_p(l,24)=v_t(i,l);
            N_p(l,7)=N_p(l,7)+v_t(i,l)*dt; %% grain radius [um]
            Cons_k(i,l)=2.102*((mean(DiffInfoD{l}(i,:,1))-xC_A_eq(i))/(xC_F_eq(i)-mean(DiffInfoD{l}(i,:,1))))^0.5871; % prefactor of Zener growth
            if t==N_p(l,6)
                v_0D(i,l)=0;
            else
                v_0D(i,l)=(Cons_k(i,l)/2)*(D_C(i)/(t-N_p(l,6)))^0.5;
            end
            v_tD(i,l)=u*v_0D(i,l);
            N_pD(l,24)=v_tD(i,l);
            N_pD(l,7)=N_pD(l,7)+v_tD(i,l)*dt;
         end
             if ll==length(Nucleated{i}(:,1))
                 for aa=1:length(A)
                      SoftPair{aa}=[];
                      HardPair{aa}=[];
                      matrixC=[];
                     for bb=1:length(A_PR{aa})
                         if A_PR{aa}(bb)~=0 && N_p(A_PR{aa}(bb),9)==1
                             bb_dis=0;
                             dis_ferrite=[];
                             GrowDis=[];
                             DiffDis=[];
                         for cc=1:length(A_PR{aa})
                             if cc~=bb && A_PR{aa}(cc)~=0 && N_p(A_PR{aa}(cc),9)==1
                                bb_dis=bb_dis+1;
                                dis_ferrite(bb_dis)=sqrt((N_p(A_PR{aa}(cc),2)-N_p(A_PR{aa}(bb),2))^2+(N_p(A_PR{aa}(cc),3)-N_p(A_PR{aa}(bb),3))^2+ ...
                                    (N_p(A_PR{aa}(cc),4)-N_p(A_PR{aa}(bb),4))^2); % distance [um]
                                if (bb_dis>1 && dis_ferrite(bb_dis)<dis_ferrite(bb_dis-1)) || bb_dis==1
                                if (mean(DiffInfo{A_PR{aa}(cc)}(i,:,4))+mean(DiffInfo{A_PR{aa}(cc)}(i,:,3))+ ...
                                        mean(DiffInfo{A_PR{aa}(bb)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,3)))<dis_ferrite(bb_dis) && N_p(A_PR{aa}(bb),25)==0
                                    N_p(A_PR{aa}(bb),25)=0; % no soft impingement
                                else if ((mean(DiffInfo{A_PR{aa}(cc)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,4)))<dis_ferrite(bb_dis) && (mean(DiffInfo{A_PR{aa}(cc)}(i,:,4))+mean(DiffInfo{A_PR{aa}(cc)}(i,:,3))+ ...
                                        mean(DiffInfo{A_PR{aa}(bb)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,3)))>=dis_ferrite(bb_dis)) || N_p(A_PR{aa}(bb),25)==1                                
                                    GrowDis(aa,A_PR{aa}(bb),A_PR{aa}(cc))=mean(DiffInfo{A_PR{aa}(bb)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,3)); % soft impingement distance [um]
                                    DiffDis(aa,A_PR{aa}(bb),A_PR{aa}(cc))=dis_ferrite(bb_dis)-mean(DiffInfo{A_PR{aa}(bb)}(i,:,4));
                                    N_p(A_PR{aa}(bb),25)=1; % soft impingement
                                    SoftPair{aa}=[SoftPair{aa};A_PR{aa}(bb) A_PR{aa}(cc) dis_ferrite(bb_dis)]; % pair of soft impingement
                                    if ((mean(DiffInfo{A_PR{aa}(cc)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,4)))>=dis_ferrite(bb_dis)) || (N_p(A_PR{aa}(bb),25)==1 && N_p(A_PR{aa}(bb),26)==1)
                                        N_p(A_PR{aa}(bb),25)=1; % soft flag
                                        N_p(A_PR{aa}(bb),26)=1; % hard impingement flag
                                        HardPair{aa}=[HardPair{aa};A_PR{aa}(bb) A_PR{aa}(cc) dis_ferrite(bb_dis)]; % pair of hard impingement
                                    end
                                    end
                                end
                                end
                             end
                         end
                         
                         if isempty(dis_ferrite) && N_p(A_PR{aa}(bb),31)~=0
                             cc=N_p(A_PR{aa}(bb),31);
                             bb_dis=bb_dis+1;
                                dis_ferrite(bb_dis)=2*(3*N_p(A_PR{aa}(bb),29)/(4*pi))^(1/3); % Local average growing d
                                if (bb_dis>1 && dis_ferrite(bb_dis)<dis_ferrite(bb_dis-1)) || bb_dis==1
                                if (mean(DiffInfo{cc}(i,:,4))+mean(DiffInfo{cc}(i,:,3))+ ...
                                        mean(DiffInfo{A_PR{aa}(bb)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,3)))<dis_ferrite(bb_dis) && N_p(A_PR{aa}(bb),25)==0
                                    N_p(A_PR{aa}(bb),25)=0; % no soft impingement
                                else if ((mean(DiffInfo{cc}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,4)))<dis_ferrite(bb_dis) && (mean(DiffInfo{cc}(i,:,4))+mean(DiffInfo{cc}(i,:,3))+ ...
                                        mean(DiffInfo{A_PR{aa}(bb)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,3)))>=dis_ferrite(bb_dis)) || N_p(A_PR{aa}(bb),25)==1
                                    
                                    GrowDis(aa,A_PR{aa}(bb),cc)=mean(DiffInfo{A_PR{aa}(bb)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,3)); % soft impingement distance [um]
                                    DiffDis(aa,A_PR{aa}(bb),cc)=N_p(A_PR{aa}(bb),32)-mean(DiffInfo{A_PR{aa}(bb)}(i,:,4)); % diffusion length does not change any more [um]
                                    N_p(A_PR{aa}(bb),25)=1; % soft impingement
                                    SoftPair{aa}=[SoftPair{aa};A_PR{aa}(bb) cc dis_ferrite(bb_dis)]; % pair of soft impingement
                                    if ((mean(DiffInfo{cc}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,4)))>=dis_ferrite(bb_dis)) || (N_p(A_PR{aa}(bb),25)==1 && N_p(A_PR{aa}(bb),26)==1)
                                        N_p(A_PR{aa}(bb),25)=1; % soft flag
                                        N_p(A_PR{aa}(bb),26)=1; % hard impingement flag
%                                         N_p(A_PR{aa}(bb),27)=mean(DiffInfo{A_PR{aa}(bb)}(i,:,4))+mean(DiffInfo{A_PR{aa}(bb)}(i,:,3)); % soft impingement distance [um]
                                        HardPair{aa}=[HardPair{aa};A_PR{aa}(bb) cc dis_ferrite(bb_dis)]; % pair of hard impingement
                                    end
                                    end
                                end
                                end
                         end
                            matrixC(aa,A_PR{aa}(bb))=mean(DiffInfo{A_PR{aa}(bb)}(i,:,1)); % matrix C [mol%]
                             if ~isempty(GrowDis) && ~isempty(DiffDis)
                                N_p(A_PR{aa}(bb),32)=min([N_p(A_PR{aa}(bb),32) min(nonzeros(GrowDis(aa,A_PR{aa}(bb),:)))]);
                                DiffInfo{A_PR{aa}(bb)}(i,:,3)=N_p(A_PR{aa}(bb),32)-mean(DiffInfo{A_PR{aa}(bb)}(i,:,4)); % diffusion length does not change any more [um]
                                if DiffInfo{A_PR{aa}(bb)}(i,:,3)<0
                                    DiffInfo{A_PR{aa}(bb)}(i,:,3)=0;
                                end
                             end
                         end
                     end
                        if ~isempty(matrixC)
                         A_P(aa,15)=mean(nonzeros(matrixC(aa,:))); % matrix C [mol%]
                         A_P(aa,16)=100*(100-A_P(aa,15))*Ux/(A_P(aa,15)+(100-A_P(aa,15))*Ux+(100-A_P(aa,15))*(1-Ux)); % Mn concentent [mol%]
                         A_P(aa,10)=100*A_P(aa,15)*M_C/(A_P(aa,15)*M_C+(100-A_P(aa,15))*Ux*M_Mn+(100-A_P(aa,15))*(1-Ux)*M_Fe); % maxtrix C [wt.%]
                         A_P(aa,11)=100*(100-A_P(aa,15))*Ux*M_Mn/(A_P(aa,15)*M_C+(100-A_P(aa,15))*Ux*M_Mn+(100-A_P(aa,15))*(1-Ux)*M_Fe); % maxtrix Mn [wt.%]
                        end
                 end
             end                    
      else
          N_p(l,7)=N_p(l,7);
      end
        ll=ll+1;
    end
    end
  
    %%%%[Vol]=impingement2(N_p,numb);%call the impingement2 function
    %[Vol Tri impinge4_flag]=impingement34(N_p,numb);%call the impingement34 function
    %Tri
    mm=1;p4=0;p3=0;delta(i)=0;volume=0;
    if ~isempty(Nucleated{i})
    % Restore the coordinates and radious of particle and its mirror
    % partice: an array of [5*numb,4]
    for j=1:length(Nucleated{i}(:,1))
        xx(j,1)=N_p(Nucleated{i}(j,1),2);
        xx(j,2)=N_p(Nucleated{i}(j,1),3);
        xx(j,3)=N_p(Nucleated{i}(j,1),4);
        xx(j,4)=N_p(Nucleated{i}(j,1),7);
        xx(length(Nucleated{i}(:,1))+j,1)=position_MXY(Nucleated{i}(j,1),1);
        xx(length(Nucleated{i}(:,1))+j,2)=position_MXY(Nucleated{i}(j,1),2);
        xx(length(Nucleated{i}(:,1))+j,3)=position_MXY(Nucleated{i}(j,1),3);
        xx(length(Nucleated{i}(:,1))+j,4)=N_p(Nucleated{i}(j,1),7);
        
        xx(length(Nucleated{i}(:,1))*2+j,1)=position_MXZ(Nucleated{i}(j,1),1);
        xx(length(Nucleated{i}(:,1))*2+j,2)=position_MXZ(Nucleated{i}(j,1),2);
        xx(length(Nucleated{i}(:,1))*2+j,3)=position_MXZ(Nucleated{i}(j,1),3);
        xx(length(Nucleated{i}(:,1))*2+j,4)=N_p(Nucleated{i}(j,1),7);
        
        xx(length(Nucleated{i}(:,1))*3+j,1)=position_MYZ(Nucleated{i}(j,1),1);
        xx(length(Nucleated{i}(:,1))*3+j,2)=position_MYZ(Nucleated{i}(j,1),2);
        xx(length(Nucleated{i}(:,1))*3+j,3)=position_MYZ(Nucleated{i}(j,1),3);
        xx(length(Nucleated{i}(:,1))*3+j,4)=N_p(Nucleated{i}(j,1),7);
        
        xx(length(Nucleated{i}(:,1))*4+j,1)=position_MO(Nucleated{i}(j,1),1);
        xx(length(Nucleated{i}(:,1))*4+j,2)=position_MO(Nucleated{i}(j,1),2);
        xx(length(Nucleated{i}(:,1))*4+j,3)=position_MO(Nucleated{i}(j,1),3);
        xx(length(Nucleated{i}(:,1))*4+j,4)=N_p(Nucleated{i}(j,1),7);
    end
    
%     [Vol Tri impinge4_flag]=impingement34_PBC(xx,length(Nucleated{i}(:,1))*5);% call the impingement34_PBC function, slower
    [Vol Tri impinge4_flag]=impingement34_PBC_modify(xx,length(Nucleated{i}(:,1))*5, ...
        [Nucleated{i}(:,13),Nucleated{i}(:,13),Nucleated{i}(:,13), ...
        Nucleated{i}(:,13),Nucleated{i}(:,13)]);% call the impingement34_PBC function, faster
    
    % Calculate the substracted volume for particle i
    for n=1:numb       
        Vol_PBC(n)=Vol(n);
        if impinge4_flag(n)==1||impinge4_flag(n+numb)==1||impinge4_flag(n+numb*2)==1||impinge4_flag(n+numb*3)==1||impinge4_flag(n+numb*4)==1
           impinge4_flag_PBC(n)=1;
        else
           impinge4_flag_PBC(n)=0;
        end
    end   
    while mm<=length(Nucleated{i}(:,1)) %calculate the mean value and the volume fraction
        m=Nucleated{i}(mm,1);
        if N_p(m,7)>0
        N_p(m,8)=(3*Vol_PBC(mm)/(4*pi))^(1/3);%calculate the effective radius
        N_p(m,13)=impinge4_flag_PBC(mm);% Transfer the impinge3_flag to the matrix N_p
        if N_p(m,13)==1
%            N_p(m,8)=((3*(1-exp(-4/3*pi*N_p(m,7)^3/N_p(m,29)))*N_p(m,29))/(4*pi))^(1/3); % correctiong of extended volume based on Avrami approach
           N_p(m,8)=((3*tanh(4/3*pi*N_p(m,7)^3/N_p(m,29))*N_p(m,29))/(4*pi))^(1/3); % correction with tanh(xe)
           if N_p(m,8)<Nucleated{i}(mm,8)
               N_p(m,8)=Nucleated{i}(mm,8)+eps/2;
           end
           N_p(m,30)=t; % restore the moment when it includes correction of extended volume
        end
        if N_p(m,8)<(N_p(m,7)-eps)
            N_p(m,25)=1;
            N_p(m,26)=1;
        end
        N_p(m,10)=4/3*pi*N_p(m,8)^3;%actual volume calculate from the effective radius
        N_p(m,11)=N_p(m,10)/(4/3*pi*N_p(m,7)^3);%actual volume/original volume,should always be <=1
        p3=p3+N_p(m,8)^3;
        p4=p4+N_p(m,8)^4;        
        volume=volume+4/3*pi*N_p(m,8)^3; %used for the impingement34 function
        end
        mm=mm+1;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update details of each austenite
     for q=1:length(A_P)
         A_P(q,6)=0;
         A_P(q,12)=0;
         A_P(q,14)=0;
        for p=1:length(N_p(:,1))
         if N_p(p,6)>=0 && N_p(p,9)==1 %%% nucleated
             for r=1:length(N_PR{p}(:,1))
                 if N_PR{p}(r,1)==q
                    A_P(q,6)=A_P(q,6)+1; % update the counts of the nucleated ferrite in each austenite
                    A_P(q,12)=A_P(q,12)+N_p(p,10)/length(N_PR{p}(:,1)); % update the ferrite volume in each austenite [um3]
                    if DiffInfo{p}(i,r,7)==1
                        A_P(q,14)=A_P(q,14)+1; % counts of soft impingement
                    end
                 end
             end
         end
        end
        A_P(q,13)=A_P(q,12)/A_P(q,5); % ferrite fraction
        if A_P(q,13)>1
            A_P(q,13)=1;
        end
     end
     AP_track{i}=A_P;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% update details of ferrite
    for s=1:length(site)
        N_p(s,15)=mean(A_P(N_PR{s}(:,1),10)); % average remote C content in austenite [wt.%]
        N_p(s,16)=(100-N_p(s,15))*Ux/(N_p(s,15)+(100-N_p(s,15))*Ux+(100-N_p(s,15))*(1-Ux))*100; % average remote Mn content in austenite [wt.%]
        N_p(s,17)=wC_F_eq(i); % average remote C content in ferrite [wt.%]
        N_p(s,18)=(100-N_p(s,17))*Ux/(N_p(s,17)+(100-N_p(s,17))*Ux+(100-N_p(s,17))*(1-Ux))*100; % average remote Mn content in ferrite [wt.%]
    end
    for Ia=1:length(Nucleated{i}(:,1)) %to exlude the nucleation position which has been contained in the ferrite and calculate the C concentration for each site
         a=Nucleated{i}(Ia,1);
        for b=1:length(site)
            if a~=b
                d(a,b)=sqrt((N_p(b,2)-N_p(a,2))^2+(N_p(b,3)-N_p(a,3))^2+(N_p(b,4)-N_p(a,4))^2);
                d_MXY(a,b)=sqrt((N_p(b,2)-position_MXY(a,1))^2+(N_p(b,3)-position_MXY(a,2))^2+(N_p(b,4)-position_MXY(a,3))^2);
                d_MXZ(a,b)=sqrt((N_p(b,2)-position_MXZ(a,1))^2+(N_p(b,3)-position_MXZ(a,2))^2+(N_p(b,4)-position_MXZ(a,3))^2);
                d_MYZ(a,b)=sqrt((N_p(b,2)-position_MYZ(a,1))^2+(N_p(b,3)-position_MYZ(a,2))^2+(N_p(b,4)-position_MYZ(a,3))^2);
                d_MO(a,b)=sqrt((N_p(b,2)-position_MO(a,1))^2+(N_p(b,3)-position_MO(a,2))^2+(N_p(b,4)-position_MO(a,3))^2);
                Rab(a,b)=N_p(a,7);
                if (d(a,b)<=N_p(a,7)&&N_p(b,9)==0)||(d_MXY(a,b)<=N_p(a,7)&&N_p(b,9)==0)|| ...
                        (d_MXZ(a,b)<=N_p(a,7)&&N_p(b,9)==0)||(d_MYZ(a,b)<=N_p(a,7)&&N_p(b,9)==0)|| ...
                        (d_MO(a,b)<=N_p(a,7)&&N_p(b,9)==0)
                   N_p(b,5)=0; % in the range of Rbcc+shield distance
                end
                for c=1:length(N_PR{a}(:,1))
                    if (d(a,b)<=(N_p(a,7)+DiffInfo{a}(i,c,3))&&N_p(b,6)==0)||(d_MXY(a,b)<=(N_p(a,7)+DiffInfo{a}(i,c,3))&&N_p(b,6)==0)|| ...
                        (d_MXZ(a,b)<=(N_p(a,7)+DiffInfo{a}(i,c,3))&&N_p(b,6)==0)||(d_MYZ(a,b)<=(N_p(a,7)+DiffInfo{a}(i,c,3))&&N_p(b,6)==0)|| ...
                        (d_MO(a,b)<=(N_p(a,7)+DiffInfo{a}(i,c,3))&&N_p(b,6)==0)
                    Lmin=min([d(a,b) d_MXY(a,b) d_MXZ(a,b) d_MYZ(a,b) d_MO(a,b)])-N_p(a,7); % minimum distance from the interface [um]
                    N_p(b,22)=DiffInfo{a}(i,c,1)+(DiffInfo{a}(i,c,2)-DiffInfo{a}(i,c,1))*(1-Lmin/DiffInfo{a}(i,c,3))^2; % carbon concentration [mol%]
                    N_p(b,23)=100*(100-N_p(b,22))*Ux*M_Mn/(N_p(b,22)*M_C+(100-N_p(b,22))*Ux*M_Mn+(100-N_p(b,22))*(1-Ux)*M_Fe); % Mn concentent [wt.%]
                    N_p(b,22)=100*N_p(b,22)*M_C/(N_p(b,22)*M_C+(100-N_p(b,22))*Ux*M_Mn+(100-N_p(b,22))*(1-Ux)*M_Fe); % carbon concentration [wt.%]
                    end
                end
            else
               N_p(b,22)=wC_F_eq(i); % C content [wt.%]
               N_p(b,23)=100*(100-N_p(b,22))*Ux/(N_p(b,22)+(100-N_p(b,22))*Ux+(100-N_p(b,22))*(1-Ux)); % Mn content [wt.%]
            end
     end
        %%%%%%%%%%%%
        N_p(a,29)=0;
        for m=1:length(N_PR{a}(:,1))
            N_p(a,29)=N_p(a,29)+A_P(N_PR{a}(m,1),5)*(1-A_P(N_PR{a}(m,1),13)); % volume that it can maximum grow [um^3]
        end
        N_p(a,29)=N_p(a,29)+N_p(a,10); % the rest surrounding austenite volume+its own volume [um^3]
     end
    
     % minimum distance from a nucleus[um]
    for b=1:length(site)
        if N_p(b,9)==0
        N_p(b,28)=min([min(nonzeros(d(:,b))) min(nonzeros(d_MXY(:,b))) min(nonzeros(d_MXZ(:,b))) min(nonzeros(d_MYZ(:,b))) min(nonzeros(d_MO(:,b)))]);
        N_p(b,12)=0; %ratio of volume fraction to equilibrium volume fraction predicted by phase diagram
        else
            N_p(b,28)=0;
            N_p(b,12)=N_p(b,10)/(Lb^3*F_eq(i)/N0);
            if N_p(b,12)>F_eq(i)
                N_p(b,12)=F_eq(i);
            end
            if ~isempty(Nucleated{i}) && length(Nucleated{i}(:,1))>=min(nonzeros(sumN))
            if N_p(b,25)==0 && N_p(b,26)==0
                minD_d=d(:,b)-Rab(:,b);
                minD_dMXY=d_MXY(:,b)-Rab(:,b);
                minD_dMXZ=d_MXZ(:,b)-Rab(:,b);
                minD_dMYZ=d_MYZ(:,b)-Rab(:,b);
                minD_dMO=d_MO(:,b)-Rab(:,b);
                minD=min([min(nonzeros(minD_d)) min(nonzeros(minD_dMXY)) min(nonzeros(minD_dMXZ)) ...
                    min(nonzeros(minD_dMYZ)) min(nonzeros(minD_dMO))]); % update the diffusion distance [um]             
            if ~isempty(find(minD_d==minD))
                N_p(b,31)=find(minD_d==minD);
            else if ~isempty(find(minD_dMXY==minD))
                    N_p(b,31)=find(minD_dMXY==minD);
                else if ~isempty(find(minD_dMXZ==minD))
                        N_p(b,31)=find(minD_dMXZ==minD);
                    else if ~isempty(find(minD_dMYZ==minD))
                            N_p(b,31)=find(minD_dMYZ==minD);
                        else if ~isempty(find(minD_dMO==minD))
                                N_p(b,31)=find(minD_dMO==minD);
                            end
                        end
                    end
                end
            end
            if N_p(b,31)~=0
                N_p(b,27)=2*((3*Lb^3/N0)/(4*pi))^(1/3)-N_p(N_p(b,31),7);
            else
                N_p(b,27)=2*((3*Lb^3/N0)/(4*pi))^(1/3);
            end
             V_Fsurrounding=0;            
             for r=1:length(N_PR{b}(:,1))
                 for rr=1:length(A_PR{N_PR{b}(r,1)})
                     if ~isempty(find(Nucleated{i}(:,1)==A_PR{N_PR{b}(r,1)}(rr))) && A_PR{N_PR{b}(r,1)}(rr)~=b
                        V_Fsurrounding=V_Fsurrounding+N_p(A_PR{N_PR{b}(r,1)}(rr),10)/length(N_PR{A_PR{N_PR{b}(r,1)}(rr)}(:,1));
                     end
                 end
             end
             if V_Fsurrounding==0 && N_p(b,31)~=0
                 N_p(b,32)=min([N_p(b,32) (3*sum(A_P(N_PR{b}(:,1),5))/(4*pi))^(1/3)-N_p(N_p(b,31),7)]);
             else
                 N_p(b,32)=(3*(sum(A_P(N_PR{b}(:,1),5))-V_Fsurrounding)/(4*pi))^(1/3);
             end
            else if N_p(b,25)==1 && N_p(b,26)==0 && Nucleated{i}(find(Nucleated{i}(:,1)==N_p(b,1)),25)==1
                    N_p(b,27)=Nucleated{i}(find(Nucleated{i}(:,1)==N_p(b,1)),27);
                    N_p(b,32)=Nucleated{i}(find(Nucleated{i}(:,1)==N_p(b,1)),32);
                else if N_p(b,26)==1
                    N_p(b,27)=Nucleated{i}(find(Nucleated{i}(:,1)==N_p(b,1)),27);
                    N_p(b,32)=Nucleated{i}(find(Nucleated{i}(:,1)==N_p(b,1)),32);
                    end
                end
                N_p(b,31)=Nucleated{i}(find(Nucleated{i}(:,1)==N_p(b,1)),31);
            end
             N_pD(b,25:32)=N_p(b,25:32);
            end
        end
    end
    
    % Update info for nucleated ferrite
    for b=1:length(Nucleated{i}(:,1))
        for a=1:length(N_p(:,1))
            if N_p(a,1)==Nucleated{i}(b,1)
               Nucleated{i}(b,5:end)=N_p(a,5:end);
               DiffInfo{Nucleated{i}(b,1)}(i,:,4)=Nucleated{i}(b,7);
            end
        end
    end
    end
    
%%%%%%%%% overall statistics
    if p3>0
        delta(i)=p4/p3; %effective average radius=<r^4>/<r^3>
        else if p3==0
        delta(i)=0;
        end
    end
    F(i)=volume/Lb^3; %volume fraction
    trans=N_p(:,8) ; % An intermediate array to install the effective radius
    R_sd(i)=std(trans(trans~=0));% Standard deviation of the effective radius
    
    if F(i)>=F_eq(i) && F(i)>0.1
%        F(i)=F(i)*(1-exp(-F(i)/F_eq(i)));  % Avrami correction
%        F(i)=F(i)*tanh(F(i)/F_eq(i));      % hyperbolic correction, E.J. Mittemeijer, F.Sommer, 2002.
       F(i)=F_eq(i)-eps/2;
    end
%%%%%%%%%%% the extended volume fraction
    if isempty(Nucleated{i})
        Fextend(i)=0;
    else
        Fextend(i)=sum(4/3*pi*Nucleated{i}(:,7).^3)/Lb^3;
    end 

%%%%% Following is to calculate magnetic configuration and ND characteristics at a constant applied field
if MagFlag==1
    Temp_progress=(T_a3-T(i))/(T_a3-Temp_eq(3,2))*100;
    if N_eff(i)>0 && T(i)<Tc
    [Mag_parameter0,kesa0,md_cos0,nd_angle0,nd_delta0,aax0]= ...
        Magnetic_calculation_modify_try(alpha0,beta0,MA_theta0,MA_phi0, ...
        N_p,xx,numb,F(i),T(i),Lb,Tc); % Call for the magnetic configuration function
    % Characteristics of magnetic structure
      for dd=1:aax0
        Mag_parameter(i,dd,:)=Mag_parameter0(dd,:);
        kesa(i,dd,:)=kesa0(dd,:);
        md_cos(i,dd,:)=md_cos0(dd,:);
        nd_angle(i,dd,:)=nd_angle0(dd,:);
        nd_delta(i,dd,:)=nd_delta0(dd,:);
      end
    end
end
           
%%%%%%%%%%%%%%%%%%%%%%%%%
    if T(i)<=Temp_eq(5,2)||F(i)>1.05*F_eq(i)
        stop_cycle=1;
    end
    i=i+1;
    t=t+dt;
end
%%%% save outputs
save('ferrite_3d_model_voronoin_PBC_ND_CNT.mat');
dlmwrite('myfile.txt',[Timer' T' N_eff'/(Lb*1e-6)^3 F' ...
    F_eq' delta' R_sd'],'delimiter',' ');

%%%%%%%%%% visulize results
figure('Name','Austenite geometry');
subplot(3,2,1);
V.plot('alpha', 0.05); % Adjust the transparancy
hold all;
plot3(A(:,1),A(:,2),A(:,3),'k+','LineWidth',2);% Plot the centroids
plot3(site(:,1),site(:,2),site(:,3),'r.','MarkerSize',20);% Plot the corners
plot3(Nucleated{i-1}(:,2),Nucleated{i-1}(:,3),Nucleated{i-1}(:,4),'b*','LineWidth',2);
grid off;
xlabel('x (micron)','FontSize',12);
ylabel('y (micron)','FontSize',12);
zlabel('z (micron)','FontSize',12);
axis([-Lb/2 Lb/2 -Lb/2 Lb/2 -Lb/2 Lb/2]);
box on;
% legend('+ centroids of austenite grains','. potential ferrite nucleation sites','* randomly selected nucleation sites');

subplot(3,2,2);
% plot(Timer',N_eff');
plot(T'-273,N_eff');
title('Nucleation at time t');
% xlabel('Time (seconds)');
xlabel('Temperature (^{o}C)');
ylabel('Numble of nuclei');

subplot(3,2,3);
% Plot the ferrite particles
for mm=1:length(Nucleated{i-1}(:,1))
    m=Nucleated{i-1}(mm,1);
    if N_p(m,7)>0
    ssphere(N_p(m,2),N_p(m,3),N_p(m,4),N_p(m,7),Lb);
    hold on;
    ssphere(position_MXY(m,1),position_MXY(m,2),position_MXY(m,3),N_p(m,7),Lb);
    hold on;
    ssphere(position_MXZ(m,1),position_MXZ(m,2),position_MXZ(m,3),N_p(m,7),Lb);
    hold on;
    ssphere(position_MYZ(m,1),position_MYZ(m,2),position_MYZ(m,3),N_p(m,7),Lb);
    hold on;
    ssphere(position_MO(m,1),position_MO(m,2),position_MO(m,3),N_p(m,7),Lb);
    hold on;
    title('3D visualization-Original');
    hold on;
    end
end
for mm=1:length(Nucleated{i-1}(:,1))
    m=Nucleated{i-1}(mm,1);
    if N_p(m,7)>0
        sphere0=[N_p(m,2),N_p(m,3),N_p(m,4),N_p(m,7)];
        [Plane0]=CutOffPlot(sphere0,Lb);
              hold on;
        sphere0=[position_MXY(m,1),position_MXY(m,2),position_MXY(m,3),N_p(m,7)];
        [Plane0]=CutOffPlot(sphere0,Lb);
             hold on;         
        sphere0=[position_MXZ(m,1),position_MXZ(m,2),position_MXZ(m,3),N_p(m,7)];
        [Plane0]=CutOffPlot(sphere0,Lb);
             hold on;
        sphere0=[position_MYZ(m,1),position_MYZ(m,2),position_MYZ(m,3),N_p(m,7)];
        [Plane0]=CutOffPlot(sphere0,Lb);
              hold on;
        sphere0=[position_MO(m,1),position_MO(m,2),position_MO(m,3),N_p(m,7)];
        [Plane0]=CutOffPlot(sphere0,Lb);
            hold on;
    end
end

% Plot the original austenite grains
[V,cells] = mpt_voronoi(A', 'bound', B); % Call for the mpt_voronoi
V.plot('alpha', 0.4); % Adjust the transparancy
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
zlabel('z (\mum)','FontSize',16);
set(get(gca,'xlabel'),'rotation',18);
set(get(gca,'ylabel'),'rotation',-25);
set(get(gca,'zlabel'),'rotation',90);
set(gca,'fontsize',14);
set(gca,'linewidth',2);
axis([-Lb/2 Lb/2 -Lb/2 Lb/2 -Lb/2 Lb/2]);
grid off;
hold off;
box on;

subplot(3,2,4);
for mm=1:length(Nucleated{i-1}(:,1))
    m=Nucleated{i-1}(mm,1);
    if N_p(m,8)>0
    ssphere(N_p(m,2),N_p(m,3),N_p(m,4),N_p(m,8),Lb);
        hold on;
    ssphere(position_MXY(m,1),position_MXY(m,2),position_MXY(m,3),N_p(m,8),Lb);
       hold on;
    ssphere(position_MXZ(m,1),position_MXZ(m,2),position_MXZ(m,3),N_p(m,8),Lb);
        hold on;
    ssphere(position_MYZ(m,1),position_MYZ(m,2),position_MYZ(m,3),N_p(m,8),Lb);
        hold on;
    ssphere(position_MO(m,1),position_MO(m,2),position_MO(m,3),N_p(m,8),Lb);
        hold on;
   % axis([0 126 0 126 0 126]);
    title('3D visualization-after substract the overlay');
    hold on;
    end
end
for mm=1:length(Nucleated{i-1}(:,1))
    m=Nucleated{i-1}(mm,1);
    if N_p(m,7)>0
        sphere0=[N_p(m,2),N_p(m,3),N_p(m,4),N_p(m,8)];
        [Plane0]=CutOffPlot(sphere0,Lb);
              hold on;
        sphere0=[position_MXY(m,1),position_MXY(m,2),position_MXY(m,3),N_p(m,8)];
        [Plane0]=CutOffPlot(sphere0,Lb);
             hold on;         
        sphere0=[position_MXZ(m,1),position_MXZ(m,2),position_MXZ(m,3),N_p(m,8)];
        [Plane0]=CutOffPlot(sphere0,Lb);
             hold on;
        sphere0=[position_MYZ(m,1),position_MYZ(m,2),position_MYZ(m,3),N_p(m,8)];
        [Plane0]=CutOffPlot(sphere0,Lb);
              hold on;
        sphere0=[position_MO(m,1),position_MO(m,2),position_MO(m,3),N_p(m,8)];
        [Plane0]=CutOffPlot(sphere0,Lb);
            hold on;
    end
end
[V,cells] = mpt_voronoi(A', 'bound', B); % Call for the mpt_voronoi
V.plot('alpha', 0.4); % Adjust the transparancy
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
zlabel('z (\mum)','FontSize',16);
set(get(gca,'xlabel'),'rotation',18);
set(get(gca,'ylabel'),'rotation',-25);
set(get(gca,'zlabel'),'rotation',90);
set(gca,'fontsize',14);
set(gca,'linewidth',2);
axis([-Lb/2 Lb/2 -Lb/2 Lb/2 -Lb/2 Lb/2]);
grid off;
hold off;
box on;

subplot(3,2,5);
% plot(Timer',delta','-bo',Timer',R_sd','r.-');
% plot(T'-273,delta','-bo',T'-273,R_sd','r.-');
plot(F',delta','-bo',F',R_sd','r.-');
legend('<r>','\sigma');
title('Average radius of ferrite grains');
% xlabel('Time (seconds)');
xlabel('Temperature (^{o}C)');
ylabel('Average radius (\mum)');

subplot(3,2,6);
% plot(Timer',F','-bo',Timer',F_eq','r');
plot(T'-273,F','-bo',T'-273,F_eq','r');
legend('modelling','equilibrium');
title('Volume fraction at time t');
% xlabel('Time (seconds)');
xlabel('Temperature (^{o}C)');
ylabel('Volume fraction of ferrite phase');

% TrackNo=fix(length(Nucleated{i-1})*rand(1,1));
TrackNo=1;
figure('Name', 'Track status for a particular ferrite');
subplot(2,3,1);
plot(Timer,DiffInfo{TrackNo}(:,1,4),'ko-');
hold all;
for pp=1:length(Timer)
    if ~isempty(Nucleated{pp}) && ~isempty(find(Nucleated{pp}(:,1)==TrackNo))
        plot(Timer(pp),Nucleated{pp}(find(Nucleated{pp}(:,1)==TrackNo),8),'ro-');
    end
end
hold off;
xlabel('Time (s)');
ylabel('Ferrite radius (\mum)');
set(gca,'fontsize',12);
set(gca,'linewidth',2);

subplot(2,3,2);
plot(Timer,DiffInfo{TrackNo}(:,1,3),'ro-');
xlabel('Time (s)');
ylabel('Diffusion length (\mum)');
set(gca,'fontsize',12);
set(gca,'linewidth',2);

subplot(2,3,3);
plot(Timer,DiffInfo{TrackNo}(:,1,7),'ro-');
xlabel('Time (s)');
ylabel('Flag of soft impingement');
set(gca,'fontsize',12);
set(gca,'linewidth',2);

subplot(2,3,4);
plot(Timer,DiffInfo{TrackNo}(:,1,1),'ro-',Timer,DiffInfo{TrackNo}(:,1,2),'b+-',Timer,DiffInfo{TrackNo}(:,1,6),'kd-');
legend('C_{matrix}','C^{i}','C^{eq}');
xlabel('Time (s)');
ylabel('C content (mol%)');
set(gca,'fontsize',12);
set(gca,'linewidth',2);

subplot(2,3,5);
plot(Timer,v_t(:,TrackNo),'ro-');
xlabel('Time (s)');
ylabel('V_{interface} (\mum/s)');
set(gca,'fontsize',12);
set(gca,'linewidth',2);

subplot(2,3,6);
plot(Timer,F,'ro-',Timer,Fextend,'b+-',Timer,tanh(Fextend),'md-',Timer,F_eq,'k-');
legend('fv','fv_{extended}','fv_{corrected}','feq');
xlabel('Time (s)');
ylabel('Ferrite fraction (\mum/s)');
set(gca,'fontsize',12);
set(gca,'linewidth',2);
for b=1:length(Timer)
    if ~isempty(Nucleated{b})
        GrowSpace(b)=Nucleated{b}(TrackNo,32);
    end
end

% plot ND data
if MagFlag==1
    figure;
    subplot(2,2,1);
    ds=setdiff(N_p(:,8),[0]);
    hh=hist(ds,0:round(max(ds)));
    hh=hh/length(ds);
    bar(0:round(max(ds)),hh);
    xlabel('Effective radius of ferrite (um)');
    ylabel('Numbers');
    title('Size distribution of ferrite');

    subplot(2,2,2);
    [ax h1 h2]=plotyy(Timer',delta',Timer',Mag_parameter(:,:,6));
    h1.Color='r';
    h1.LineStyle='-';
    h1.Marker='o';
    h2.Color='k';
    h2.LineStyle='-.';
    h2.Marker='d';
    hold on;
    plot(Timer',nd_delta,'b*--');
    hold off;
    title('Average radius of ferrite grains');
    legend('modelling','ND derivation');
    xlabel('Time (min)');
    ylabel('Average radius (um)');
    axes(ax(2));
    hold on;
    legend('Reduced magnetization');
    ylabel('Reduced magnetization m');
    hold off;

    subplot(2,2,3);
    [ax h1 h2]=plotyy(F',delta',F',Mag_parameter(:,:,6));
    h1.Color='r';
    h1.LineStyle='-';
    h1.Marker='o';
    h2.Color='k';
    h2.LineStyle='-.';
    h2.Marker='d';
    hold on;
    plot(F',nd_delta,'b*--');
    hold off;
    title('Average radius of ferrite grains');
    legend('modelling','ND derivation');
    xlabel('Volume fraction');
    ylabel('Average radius (um)');
    axes(ax(2));
    hold on;
    legend('Reduced magnetization');
    ylabel('Reduced magnetization m');
    hold off;

    subplot(2,2,4);
    plot(Timer',nd_angle,'ro');
    xlabel('Time (min)');
    ylabel('Depolarization rotation (rad)');
end



