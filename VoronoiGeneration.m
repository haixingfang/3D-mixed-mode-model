% This is a function to generate a cubic box containing the Voronoi cell representing the austenite grain 
% and the grain corner marked as potential nucleation sites for ferrite

% This function requires the installation of mpt3.m which can be downloaded
% at http://control.ee.ethz.ch/~mpt/3/Main/Installation;
% To add this toolbox, just type the commands: 
% tbxmanager restorepath;mpt_init;

% input parameters
clear all;
% tbxmanager restorepath;
% mpt_init;
Npot=50;     %number of potential nucleation sites [-]
f_N=1;       %scaling factor
N0=Npot*f_N; %scaled number of potential nucleation sites [-]
Lb=70;       %box size half from 126 um to 63 um, means number of austenite decreased from 250 to 31 [um]
dmin=12;     % The minimum distance required between the centroids of austenite grains
rou=42/Lb^3; % Number density of austenite grain [um-3]
enlarge_Lb=Lb+2*Lb/250^(1/3);
enlarge_N=round(rou*(Lb+2*Lb/250^(1/3))^3); % Number of austenite grains in enlargement box (enlarged by 2*average neibouring distance) that have the same number density
minus_edge=Lb/250^(1/3);

%Generate positions potential nucleation sites and plot the 3d voronoi
%diagram
%A=Lb*rand([250,3]);% Centroid coordinates of the 250 austenite grains
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
dis=pdist(A,'euclidean'); % Pairwise distance of austenite grain center

% Following is only used to validate the minimum distance follows the
% criterion
% for i=1:length(A)
%     for j=1:length(A)
%         if i~=j
%             dis(i,j)=sqrt((A(i,1)-A(j,1))^2+(A(i,2)-A(j,2))^2+(A(i,3)-A(j,3))^2);
%         else
%             dis(i,j)=inf;
%         end
%     end
% end
%plot3(A(:,1),A(:,2),A(:,3),'r.')
%B=Polyhedron([0 0 0;enlarge_Lb 0 0;enlarge_Lb enlarge_Lb 0;0 enlarge_Lb 0; 0 0 enlarge_Lb;enlarge_Lb 0 enlarge_Lb;enlarge_Lb enlarge_Lb enlarge_Lb;0 enlarge_Lb enlarge_Lb]);% Boundary vertices
B=Polyhedron([-minus_edge-Lb/2 -minus_edge-Lb/2 -minus_edge-Lb/2; ...
    enlarge_Lb-minus_edge-Lb/2 -minus_edge-Lb/2 -minus_edge-Lb/2; ...
    enlarge_Lb-minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2 -minus_edge-Lb/2; ...
    -minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2 -minus_edge-Lb/2; ...
    -minus_edge-Lb/2 -minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2; ...
    enlarge_Lb-minus_edge-Lb/2 -minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2; ...
    enlarge_Lb-minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2; ...
    -minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2]);% Boundary vertices
%B=Polyhedron([0 0 0;126 0 0;126 126 0;0 126 0; 0 0 126;126 0 126;126 126 126;0 126 126]);
figure;
[V,cells] = mpt_voronoi(A', 'bound', B); % Call for the mpt_voronoi
V.plot('alpha', 0.13); % Adjust the transparancy
hold on;
plot3(A(:,1),A(:,2),A(:,3),'k+','LineWidth',2);% Plot the centroids
vertices = V.Set.forEach(@(e) e.V, 'UniformOutput', false);% find the vertices for each cell
new_vertices = cat(1, vertices{:}); % Combine the vertices into one matrix
potential0=unique(round(new_vertices*1e6), 'rows')/1e6; % Get rid of the numerical noise 
to_remove = [-Lb/2 -Lb/2 -Lb/2;Lb/2 -Lb/2 -Lb/2;Lb/2 Lb/2 -Lb/2;-Lb/2 Lb/2 -Lb/2; ...
    -Lb/2 -Lb/2 Lb/2;Lb/2 -Lb/2 Lb/2;Lb/2 Lb/2 Lb/2;-Lb/2 Lb/2 Lb/2];% Remove the corners of the Lb*Lb*Lb cubic
potential=setdiff(potential0, to_remove, 'rows');
%potential=potential0;
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
  %if all(potential(i,:))&&all(potential(i,:)-Lb-1)
  if all(potential(i,:)-Lb/2-1)&&all(potential(i,:)+Lb/2+1)
     eval(['site',num2str(i),'=','potential(i,:)']); % Convert the number to string
     eval(['site=[site;site',num2str(i),'];']);  % Combine site1, site 2,...
  end
end
plot3(site(:,1),site(:,2),site(:,3),'r.','MarkerSize',20);
p1=randperm(length(site));
position=site(p1(1:N0),:);% Randomly select N0 potential nucleation site for ferrite
plot3(position(:,1),position(:,2),position(:,3),'b*','LineWidth',2);
grid off;
xlabel('x (micron)','FontSize',12);
ylabel('y (micron)','FontSize',12);
zlabel('z (micron)','FontSize',12);
%axis([-minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2 -minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2 -minus_edge-Lb/2 enlarge_Lb-minus_edge-Lb/2]);
axis([-Lb/2 Lb/2 -Lb/2 Lb/2 -Lb/2 Lb/2]);
box on;
legend('+ centroids of austenite grains','. potential ferrite nucleation sites');
hold off;

