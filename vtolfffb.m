function dx = vtolfffb(t,x,f,U,k,K1)
lamda=1;eplsion=1;M=1;J=1;g=1;
u=U(1,:);u1=U(2,:);u2=U(3,:);
% yd=k(1,:);ydd=k(2,:);
u1=interp1(u,u1,t);u2=interp1(u,u2,t);
% yd=interp1(u,yd,t);ydd=interp1(u,ydd,t);
x1=k(1,:);x2=k(2,:);x3=k(3,:);x4=k(4,:);x5=k(5,:);x6=k(6,:);
x1=interp1(u,x1,t);x2=interp1(u,x2,t);x3=interp1(u,x3,t);x4=interp1(u,x4,t);x5=interp1(u,x5,t);x6=interp1(u,x6,t);
% xref=[yd;ydd;yd;ydd;0;0];
xref=[x1;x2;x3;x4;x5;x6];
u1 = u1+K1(1,:)*xref;
u2 = u2+K1(2,:)*xref;
dx=zeros(6,1);
dx(1)=x(2);
dx(2)=-sin(x(5))/M*(u1-K1(1)*x(1)-K1(3)*x(2)-K1(5)*x(3)-K1(7)*x(4)-K1(9)*x(5)-K1(11)*x(6))+eplsion*cos(x(5))/M*(u2-K1(2)*x(1)-K1(4)*x(2)-K1(6)*x(3)-K1(8)*x(4)-K1(10)*x(5)-K1(12)*x(6));
dx(3)=x(4);
dx(4)=-g+cos(x(5))/M*(u1-K1(1)*x(1)-K1(3)*x(2)-K1(5)*x(3)-K1(7)*x(4)-K1(9)*x(5)-K1(11)*x(6))+eplsion*sin(x(5))/M*(u2-K1(2)*x(1)-K1(4)*x(2)-K1(6)*x(3)-K1(8)*x(4)-K1(10)*x(5)-K1(12)*x(6));
dx(5)=x(6);
dx(6)=lamda/J*(u2-K1(2)*x(1)-K1(4)*x(2)-K1(6)*x(3)-K1(8)*x(4)-K1(10)*x(5)-K1(12)*x(6));
end