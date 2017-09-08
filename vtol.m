function dx= vtol(t,x,f,U)
lamda=1;eplsion=1;M=1;J=1;g=1;
u=U(1,:);u1=U(2,:);u2=U(3,:);
u1=interp1(u,u1,t);u2=interp1(u,u2,t);
% u1=U(1,:);u2=U(2,:);
dx=zeros(6,1);
dx(1)=x(2);
dx(2)=-sin(x(5))/M*u1+eplsion*cos(x(5))/M*u2;
dx(3)=x(4);
dx(4)=-g+cos(x(5))/M*u1+eplsion*sin(x(5))/M*u2;
dx(5)=x(6);
dx(6)=lamda/J*u2;
end