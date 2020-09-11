function ff=magnetization_vector(N,easy_axis,Hl,fai_v,theta_v)
% This function is to calculation three components of the magnetization unit vector which shares the same plane with easy axis and Hk
% magnetization_vector=aa*easy_axis+bb*Hk

Hll=sqrt(Hl(1)^2+Hl(2)^2+Hl(3)^2);
Hl(1)=Hl(1)/Hll;
Hl(2)=Hl(2)/Hll;
Hl(3)=Hl(3)/Hll;% Convert Hl-vector to a unit vector

ff(1)=(N(1)*Hl(1)+N(2)*easy_axis(1))^2+(N(1)*Hl(2)+N(2)*easy_axis(2))^2+(N(1)*Hl(3)+N(2)*easy_axis(3))^2-1;% length of unit vector = 1
ff(2)=((N(1)*Hl(1)+N(2)*easy_axis(1))*Hl(1)+(N(1)*Hl(2)+N(2)*easy_axis(2))*Hl(2)+(N(1)*Hl(3)+N(2)*easy_axis(3))*Hl(3))-cos(fai_v); % angle between m and Hl
ff(3)=((N(1)*Hl(1)+N(2)*easy_axis(1))*easy_axis(1)+(N(1)*Hl(2)+N(2)*easy_axis(2))*easy_axis(2)+(N(1)*Hl(3)+N(2)*easy_axis(3))*easy_axis(3))-cos(theta_v); % angle between m and easy axis
end

