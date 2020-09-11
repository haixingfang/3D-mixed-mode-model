%This function return the volume of sphere i in which the impingement of only two spheres is considered

function [V]=impingement2(x,N0)%funtion impingement,x is the N0*7 array contains all the information,N0 is numbers of the lines
%N_p=zeros(length(x(:,1)),length(x(:,2))); %check assigned matrix size
%hatch(1)=0;
for i=1:N0
    V(i)=4/3*pi*x(i,7)^3;
    for j=1:N0
        if i~=j
            d(i,j)=((x(i,2)-x(j,2))^2+(x(i,3)-x(j,3))^2+(x(i,4)-x(j,4))^2)^0.5;%distance between grain 1 and 2
            r(i,j)=x(i,7)+x(j,7);%sum of the radius of grain 1 and 2
          if x(i,7)~=0
            if d(i,j)<r(i,j) %impingement
                cos_beta=(x(i,7)^2+d(i,j)^2-x(j,7)^2)/(2*x(i,7)*d(i,j));
                segment(i,j)=1/3*pi*x(i,7)^3*((1-cos_beta)^2)*(2+cos_beta); %overlay volume of i due to j             
             else %no impingement
                segment(i,j)=0; 
            end
          else segment(i,j)=0;
          end
        else
            segment(i,j)=0;
        end
    V(i)=V(i)-abs(segment(i,j));
       if V(i)<0
            V(i)=0;
       else V(i)=V(i);
       end
    end
    
end
