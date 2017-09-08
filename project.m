clc
clear all
close all
eplison=0.5;
A=[0 1 0 0 0 0;
    0 0 0 0 -1 0;
    0 0 0 1 0 0;
    0 0 0 0 0 0 ;
    0 0 0 0 0 1;
     0 0 0 0 0 0];
 C=[1 0 0 0 0 0;
     0 0 1 0 0 0];
% C=[1 0 0 0 0 0];
 D=0;
B=[0 0;0 -eplison;0 0;-1 0;0 0;0 1];
R=1/1000*eye(2);
Q=C'*C;
T = 2;
h=0.001;
tspan=100;
t=0:h:tspan;
yd=zeros(1,tspan/h);
for j=1:(tspan/h+1)
    if j<(T/h)
yd(:,j)=-1/(2*pi)*sin(2*pi*(j*h)/T)+1*j*h/T;
    else
        yd(:,j)=1;
    end
end
[k,a,b]=lqr(A,B,Q,R);
Acl=A-B*k;
sysclose=ss(Acl,B,C,D);
go=dcgain(sysclose);
x0=[1 0 0 0 0 0]';
r=[yd/go(1);0*yd/go(2)];

[y,t,x]=lsim(sysclose,r,t,x0);
figure(1)
plot(t,y(:,1),t,yd)