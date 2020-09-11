function [Symmetric_MXY Symmetric_MXZ Symmetric_MYZ Symmetric_MO]=SymmetricPoints(position,Lb,dmin)
% This function is to find four symmetric points of interest for original
% point position: N0*3 array
% Lb is the box size
% Returns the coordinates of the mirror points:  
% model

% Find mirror points
for i=1:length(position(:,1))
%      Symmetric_MXY(i,1)=position(i,1);
%      Symmetric_MXY(i,2)=position(i,2);
%      Symmetric_MXY(i,3)=-position(i,3);
% 
%      Symmetric_MXZ(i,1)=position(i,1);
%      Symmetric_MXZ(i,2)=-position(i,2);
%      Symmetric_MXZ(i,3)=position(i,3);
% 
%      Symmetric_MYZ(i,1)=-position(i,1);
%      Symmetric_MYZ(i,2)=position(i,2);
%      Symmetric_MYZ(i,3)=position(i,3);
% 
%      Symmetric_MO(i,1)=-position(i,1);
%      Symmetric_MO(i,2)=-position(i,2);
%      Symmetric_MO(i,3)=-position(i,3);
     
     if abs(position(i,1))<=dmin/2 && abs(position(i,2))<=dmin/2 && abs(position(i,3))<=dmin/2
         Symmetric_MXY=[];
         Symmetric_MXZ=[];
         Symmetric_MYZ=[]
         Symmetric_MO(i,1)=-position(i,1);
         Symmetric_MO(i,2)=-position(i,2);
         Symmetric_MO(i,3)=-position(i,3);
     else if abs(position(i,1))<=dmin/2 && abs(position(i,2))<=dmin/2 && abs(position(i,3))>dmin/2
         Symmetric_MXY(i,1)=position(i,1);
         Symmetric_MXY(i,2)=position(i,2);
         Symmetric_MXY(i,3)=-position(i,3);
         Symmetric_MXZ=[];
         Symmetric_MYZ=[];
         Symmetric_MO(i,1)=-position(i,1);
         Symmetric_MO(i,2)=-position(i,2);
         Symmetric_MO(i,3)=-position(i,3);
         else if abs(position(i,1))<=dmin/2 && abs(position(i,2))>dmin/2 && abs(position(i,3))<=dmin/2
         Symmetric_MXY=[];
         Symmetric_MXZ(i,1)=position(i,1);
         Symmetric_MXZ(i,2)=-position(i,2);
         Symmetric_MXZ(i,3)=position(i,3);
         Symmetric_MYZ=[];
         Symmetric_MO(i,1)=-position(i,1);
         Symmetric_MO(i,2)=-position(i,2);
         Symmetric_MO(i,3)=-position(i,3);
             else if abs(position(i,1))>dmin/2 && abs(position(i,2))<=dmin/2 && abs(position(i,3))<=dmin/2
         Symmetric_MXY=[];
         Symmetric_MXZ=[];
         Symmetric_MYZ(i,1)=-position(i,1);
         Symmetric_MYZ(i,2)=position(i,2);
         Symmetric_MYZ(i,3)=position(i,3);

         Symmetric_MO(i,1)=-position(i,1);
         Symmetric_MO(i,2)=-position(i,2);
         Symmetric_MO(i,3)=-position(i,3);
         else if abs(position(i,1))>dmin/2 && abs(position(i,2))>dmin/2 && abs(position(i,3))<=dmin/2
         Symmetric_MXY=[];
         Symmetric_MXZ(i,1)=position(i,1);
         Symmetric_MXZ(i,2)=-position(i,2);
         Symmetric_MXZ(i,3)=position(i,3);

         Symmetric_MYZ(i,1)=-position(i,1);
         Symmetric_MYZ(i,2)=position(i,2);
         Symmetric_MYZ(i,3)=position(i,3);

         Symmetric_MO(i,1)=-position(i,1);
         Symmetric_MO(i,2)=-position(i,2);
         Symmetric_MO(i,3)=-position(i,3);
                  else if abs(position(i,1))>dmin/2 && abs(position(i,2))<=dmin/2 && abs(position(i,3))>dmin/2
         Symmetric_MXY(i,1)=position(i,1);
         Symmetric_MXY(i,2)=position(i,2);
         Symmetric_MXY(i,3)=-position(i,3);

         Symmetric_MXZ=[];
         
         Symmetric_MYZ(i,1)=-position(i,1);
         Symmetric_MYZ(i,2)=position(i,2);
         Symmetric_MYZ(i,3)=position(i,3);

         Symmetric_MO(i,1)=-position(i,1);
         Symmetric_MO(i,2)=-position(i,2);
         Symmetric_MO(i,3)=-position(i,3);
                           else if abs(position(i,1))<=dmin/2 && abs(position(i,2))>dmin/2 && abs(position(i,3))>dmin/2
         Symmetric_MXY(i,1)=position(i,1);
         Symmetric_MXY(i,2)=position(i,2);
         Symmetric_MXY(i,3)=-position(i,3);
         
         Symmetric_MXZ(i,1)=position(i,1);
         Symmetric_MXZ(i,2)=-position(i,2);
         Symmetric_MXZ(i,3)=position(i,3);
         
         Symmetric_MYZ=[];
         
         Symmetric_MO(i,1)=-position(i,1);
         Symmetric_MO(i,2)=-position(i,2);
         Symmetric_MO(i,3)=-position(i,3);
         else if abs(position(i,1))>dmin/2 && abs(position(i,2))>dmin/2 && abs(position(i,3))>dmin/2         
                 Symmetric_MXY(i,1)=position(i,1);
                 Symmetric_MXY(i,2)=position(i,2);
                 Symmetric_MXY(i,3)=-position(i,3);

                 Symmetric_MXZ(i,1)=position(i,1);
                 Symmetric_MXZ(i,2)=-position(i,2);
                 Symmetric_MXZ(i,3)=position(i,3);

                 Symmetric_MYZ(i,1)=-position(i,1);
                 Symmetric_MYZ(i,2)=position(i,2);
                 Symmetric_MYZ(i,3)=position(i,3);

                 Symmetric_MO(i,1)=-position(i,1);
                 Symmetric_MO(i,2)=-position(i,2);
                 Symmetric_MO(i,3)=-position(i,3);
             end
                               end
                      end
             end
                 end
             end
         end
     end
end
end

    