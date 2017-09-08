clc
clear all
close all
% Step 1 %%%%%%%%%%%%%%%%% Define the desired output %%%%%
T = 1; % T is the time period in seconds
om = 2*pi/T; delt =0.001;
lamda=1;eplsion=0.1;M=1;J=1;g=9.81;
% the step size for computation
Ymax = 1; Ay = Ymax*2*pi/(T^2);
%
t = 0:delt:T;
yd = (Ay/om)*t -(Ay/om/om)*sin(om*t);
ydd = (Ay/om)*(1 - cos(om*t));
yd2d = Ay*sin(om*t);

f=2000;
a=2*pi*f;
Gfilter=tf([0 a],[1 a]);

z0=yd;
z1=lsim(Gfilter,z0,t);
yd=lsim(Gfilter,z1,t);
ydd=a*(z1-yd);
yd2d=a^2*(z0'-2*z1+yd);
% add zeros/constants
Npre =5; Npost =Npre;
tpre = 0:delt:(Npre*T-delt);
tpost = 0:delt:(Npost*T-delt);
%
t = 0:delt:(1+Npre+Npost)*T;
yd = [zeros(size(tpre)) yd' Ymax*ones(size(tpost))];
ydd = [zeros(size(tpre)) ydd' zeros(size(tpost))];
yd2d = [zeros(size(tpre)) yd2d' zeros(size(tpost))];
plot(t,yd,t,ydd,t,yd2d)
%%
iteration=16; % number of iterations
% Initializing
eta2_total=[]; % store unstable solution at each iteration
eta1_total=[]; % store stable solution at each iteration
etau = zeros(size(yd));
% *****************************************
% Starting the loop for the zero dynamics
figure
plot(t,yd,t,ydd,t,yd2d)

%%%calculate approxiated eta
gama=sqrt(lamda*M*g/J/eplsion);
eta=etau;
Au=gama; As=-gama; Bu=1/(2*gama); Bs=-1/(2*gama); C=1; D=0;

% Aeta=[0 1;(gama)^2 0];Beta=eye(2);
% [m,~]=eig(Aeta);Astar=inv(m)*Aeta*m;Bstar=inv(m)*Beta;
% Au=Astar(1,1);As=Astar(2,2);Bu=Bstar(1,:);Bs=Bstar(2,:);

zd2d=zeros(size(yd2d));
zd2d=yd2d;
for i=1:iteration
pn=(lamda*M*yd2d.*cos(eta)+lamda*(M*zd2d+M*g).*sin(eta)-lamda*M*g*eta)/(J*eplsion);
u_tf = -inv(Au)*Bu*pn(end);
[Zu,zu]=lsim(-Au,-Bu,C,D,fliplr(pn),t,u_tf);
Zu = fliplr(Zu');
u_to = -inv(As)*Bs*pn(1);
[Zs,zs]=lsim(As,Bs,C,D,pn,t,u_to);
Zs = Zs';
eta=inv([gama -1;gama 1])*2*gama*[Zs;Zu];
eta1_total=[eta1_total; eta(1,:)];
eta2_total=[eta2_total; eta(2,:)];
eta=eta(1,:);
end
figure(2)
subplot(211), plot(t,eta1_total')
xlabel('time')
ylabel('\eta_1')
subplot(212), plot(t,eta2_total')
xlabel('time')
ylabel('\eta_2')
%%%%%simulate
eta1=eta1_total(iteration,:);eta2=eta2_total(iteration,:);
uff=1/eplsion*[-eplsion*sin(eta1)*M.*yd2d+eplsion*cos(eta1).*(M*zd2d+M*g);cos(eta1)*M.*yd2d+sin(eta1).*(M*zd2d+M*g)];
figure(3)
plot(t,uff)
xlabel('time')
ylabel('inverse input')
legend('u_1ff,','u_2ff')

U=[t;[uff(1,:);uff(2,:)]];
x0=[0 0 0 0 0 0]';
[t,x]=ode45('vtol',t,x0,[],U);
figure(4)
plot(t,yd,t,x(:,1),t,x(:,3))
xlabel('time')
ylabel('output')
legend('desired trajectory','y1','y2')

%feedback controller
A=[0 1 0 0 0 0;
    0 0 0 0 -g 0;
    0 0 0 1 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 1;
    0 0 0 0 0 0];
B=[0 0;0 eplsion/M;0 0;1/M 0;0 0;0 lamda/J];
C=[1 0 0 0 0 0;0 0 1 0 0 0];
D=0;
R=1/1000*eye(2);Q=C'*C;
x0=[0 0 0 0 0 0]';
[K1,S,e]=lqr(A,B,Q,R);
A1=A-B*K1;
sys_close=ss(A1,B,C,D);
go= dcgain(sys_close);%calculate go
% r=[yd/go(2);-yd/go(2)];
r=inv(go)*[yd;yd];
[y1,t,x1]=lsim(sys_close,r,t,x0);%simulate 
figure(5)
plot(t,y1(:,1),t,y1(:,2),t,yd)
title('desired trajectory vs fb simulation,T=5')
xlabel('t')
ylabel('y')
legend('y1','y2','desired trajectory')
max(abs(y1(:,1)'-yd));
max(abs(y1(:,2)'-yd));

%add feedback to feedforward vs feedback alone in ODE45
U=[t';[uff(1,:);uff(2,:)]];
x0=[0;0;0;0;0;0];
k=[yd;ydd;yd;ydd;eta1;eta2];

[t3,x3]=ode45('vtolfffb',t,x0,[],U,k,K1);
figure(6)
plot(t,yd,t,x3(:,1),t,x3(:,3))
xlabel('time')
ylabel('output')
legend('desired trajectory','y1_{ff+fb}','y2_{ff+fb}')

figure(7)
plot(t,x3(:,1),t,x3(:,3),t,yd,t,y1(:,1),'--',t,y1(:,2),'--')
title('desired trajectory vs fb vs fb+ff T=5')
xlabel('t')
ylabel('y')
legend('y1_{fffb}','y2_{fffb}','desired trajectory','y1_{fb}','y2_{fb}')

ufb=-K1*(x3'-k);
figure(8)
plot(t,ufb(1,:),t,ufb(2,:),t,uff(1,:),'--',t,uff(2,:),'--')
title('ff+fb input vs fb input T=5')
xlabel('t')
ylabel('u')
legend('u1_{fb}','u2_{fb}','u1_{ff}','u2_{ff}')