clc
clear all
close all
% Step 1 %%%%%%%%%%%%%%%%% Define the desired output %%%%%
T = 10; % T is the time period in seconds
om = 2*pi/T; delt =0.001;
lamda=1;eplsion=1;M=1;J=1;g=1;
% the step size for computation
Ymax = 1; Ay = Ymax*2*pi/(T^2);
%
t = 0:delt:T;
yd = (Ay/om)*t -(Ay/om/om)*sin(om*t);zd=cos(om*t);
ydd = (Ay/om)*(1 - cos(om*t));zdd=-0.5*om*sin(om*t);
yd2d = Ay*sin(om*t);zd2d=-om*om*cos(om*t);
plot(yd,zd)

% add zeros/constants
Npre =4; Npost =Npre;
tpre = 0:delt:(Npre*T-delt);
tpost = 0:delt:(Npost*T-delt);
%
t = 0:delt:(1+Npre+Npost)*T;
yd = [zeros(size(tpre)) yd ones(size(tpost))];
ydd = [zeros(size(tpre)) ydd zeros(size(tpost))];
yd2d = [zeros(size(tpre)) yd2d zeros(size(tpost))];
%
zd = [ones(size(tpre)) zd Ymax*ones(size(tpost))];
zdd = [zeros(size(tpre)) zdd zeros(size(tpost))];
zd2d = [zeros(size(tpre)) zd2d zeros(size(tpost))];
plot(t,yd,t,ydd,t,yd2d)

iteration=16; % number of iterations
% Initializing
eta2_total=[]; eta2_hat=[];% store unstable solution at each iteration
eta1_total=[];eta1_hat=[]; % store stable solution at each iteration
etau = zeros(size(yd));


iteration=16; % number of iterations
% Initializing
eta2_total=[]; eta2_hat=[];% store unstable solution at each iteration
eta1_total=[];eta1_hat=[]; % store stable solution at each iteration
etau = zeros(size(yd));


%%%calculate approxiated eta
gama=sqrt(lamda*M*g/J/eplsion);
eta=etau;
Au=gama; As=-gama; Bu=1/(2*gama); Bs=-1/(2*gama); C=1; D=0;

%calcualte eta_stable for later use
for i=1:iteration
pn=(lamda*M*yd2d.*cos(eta)+lamda*(M*zd2d+M*g).*sin(eta)-lamda*M*g*eta)/(J*eplsion);
u_tf = -inv(Au)*Bu*pn(end);
[Zu,zu1]=lsim(-Au,-Bu,C,D,fliplr(pn),t,u_tf);
Zu = fliplr(Zu');
u_to = -inv(As)*Bs*pn(1);
[Zs,zs1]=lsim(As,Bs,C,D,pn,t,u_to);
Zs = Zs';
eta=inv([gama -1;gama 1])*2*gama*[Zs;Zu];
eta1_total=[eta1_total; eta(1,:)];
eta2_total=[eta2_total; eta(2,:)];
eta=eta(1,:);
end
%fisrt Tp*m part of smulation
Tp=0.8;m=112;it=8;
eta_0=[exp(-gama*(t(1:(Tp/delt))))*eta(Tp/delt+1);zeros(1,Tp/delt)];
for j=1:m
for i=1:it
yd2d_pre=yd2d((Tp/delt*(j-1)+1):(Tp/delt*j));
zd2d_pre=zd2d((Tp/delt*(j-1)+1):(Tp/delt*j));
t_pre=t(1:(Tp/delt));

E=exp(-gama*(t_pre))*eta(Tp/delt*(j-1)+1);
% eta_0=[E;zeros(1,length(E))];
eta_hat_pre=[1 0]*inv([gama -1;gama 1])*2*gama*eta_0;
pn=(lamda*M*yd2d_pre.*cos(eta_hat_pre)+lamda*(M*zd2d_pre+M*g).*sin(eta_hat_pre)-lamda*M*g*eta_hat_pre)/(J*eplsion);

u_tf = -inv(Au)*Bu*pn(end);
[Zu,zu2]=lsim(-Au,-Bu,C,D,fliplr(pn),t_pre,u_tf);
Zu = fliplr(Zu');

u_to = -inv(As)*Bs*pn(1);
[Zs,zs2]=lsim(As,Bs,C,D,pn,t_pre,u_to);
Zs = Zs'+E;
eta_0=[Zs;Zu];
end
eta_f=inv([gama -1;gama 1])*2*gama*[Zs;Zu];
eta1_hat=[eta1_hat eta_f(1,:)];
eta2_hat=[eta2_hat eta_f(2,:)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rest of t-Tp*m part of simulation
eta_0=eta_0(:,(Tp/delt+1-(length(t)-Tp*m/delt)):end);
for i=1:it
yd2d_pre=yd2d(Tp*m/delt+1:length(t));
zd2d_pre=zd2d(Tp*m/delt+1:length(t));
t_pre=t(1:(length(t)-Tp*m/delt));

E=exp(-gama*(t_pre))*eta(Tp*m/delt+1);
eta_hat_pre=[1 0]*inv([gama -1;gama 1])*2*gama*eta_0;
pn=(lamda*M*yd2d_pre.*cos(eta_hat_pre)+lamda*(M*zd2d_pre+M*g).*sin(eta_hat_pre)-lamda*M*g*eta_hat_pre)/(J*eplsion);

u_tf = -inv(Au)*Bu*pn(end);
[Zu,zu]=lsim(-Au,-Bu,C,D,fliplr(pn),t_pre,u_tf);
Zu = fliplr(Zu');

u_to = -inv(As)*Bs*pn(1);
u_to = Zs(end);
[Zs,zs]=lsim(As,Bs,C,D,pn,t_pre,u_to);
Zs = Zs'+E;
eta_0=[Zs;Zu];
end
eta_f=inv([gama -1;gama 1])*2*gama*[Zs;Zu];
eta1_hat=[eta1_hat eta_f(1,:)];
eta2_hat=[eta2_hat eta_f(2,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%simulation
eta1=eta1_hat;eta2=eta2_hat;
uff=1/eplsion*[-eplsion*sin(eta1)*M.*yd2d+eplsion*cos(eta1).*(M*zd2d+M*g);cos(eta1)*M.*yd2d+sin(eta1).*(M*zd2d+M*g)];
figure(3)
plot(t,uff)
xlabel('time')
ylabel('inverse input')
legend('u_1ff,','u_2ff')

U=[t;[uff(1,:);uff(2,:)]];
x0=[0 0 1 0 0 0]';
[t,x]=ode45('vtol',t,x0,[],U);

figure(4)
subplot(211)
plot(t,yd,t,x(:,1))
xlabel('time')
ylabel('output')
legend('desired trajectory x','tracking trajectory x_{ff}')
subplot(212)
plot(t,zd,t,x(:,3))
xlabel('time')
ylabel('output')
legend('desired trajectory z','tracking trajectory z_{ff}')

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
R=1/100000*eye(2);Q=C'*C;
[K1,S,e]=lqr(A,B,Q,R);
A1=A-B*K1;
sys_close=ss(A1,B,C,D);
go= dcgain(sys_close);%calculate go
r=inv(go)*[yd;zd];
[y1,t,x1]=lsim(sys_close,r,t,x0);%simulate 

figure(5)
subplot(211)
plot(t,y1(:,1),t,yd)
title('desired trajectory vs fb simulation,T=5')
xlabel('t')
ylabel('x')
legend('tracking trajectory x_{fb}','desired trajectory')
subplot(212)
plot(t,y1(:,2),t,zd)
title('desired trajectory vs fb simulation,T=5')
xlabel('t')
ylabel('z')
legend('tracking trajectory z_{fb}','desired trajectory')


%add feedback to feedforward vs feedback alone in ODE45
U=[t';[uff(1,:);uff(2,:)]];
k=[yd;ydd;zd;zdd;eta1;eta2];

[t3,x3]=ode45('vtolfffb',t,x0,[],U,k,K1);
figure(6)
subplot(211)
plot(t,yd,t,x3(:,1))
xlabel('time')
ylabel('output x')
legend('desired trajectory x','x_{ff+fb}')
subplot(212)
plot(t,zd,t,x3(:,3))
xlabel('time')
ylabel('output z')
legend('desired trajectory z','z_{ff+fb}')

figure(7)
subplot(211)
plot(t,x3(:,1),t,yd,t,y1(:,1),'--')
title('desired trajectory vs fb vs fb+ff T=5')
xlabel('t')
ylabel('x')
legend('x_{fffb}','desired trajectory','x_{fb}')
subplot(212)
plot(t,x3(:,3),t,zd,t,y1(:,2),'--')
title('desired trajectory vs fb vs fb+ff T=5')
xlabel('t')
ylabel('z')
legend('z_{fffb}','desired trajectory','z_{fb}')


ufb=-K1*(x3'-k);
figure(8)
plot(t,ufb(1,:),t,ufb(2,:),t,uff(1,:),'--',t,uff(2,:),'--')
title('ff+fb input vs fb input T=5')
xlabel('t')
ylabel('u')
legend('u1_{fb}','u2_{fb}','u1_{ff}','u2_{ff}')

figure(14)
plot(yd,zd)
hold on
plot(x3(:,1),x3(:,3),'--')
hold on
plot(x1(:,1),x1(:,3))
legend('desired','ff+fb','fb')
title('tracking trajectory vs desired trajectory, Tp=0.8s')

figure(15)
subplot(211)
plot(t,yd,t,x1(:,1),t,x3(:,1))
xlabel('t')
ylabel('x')
legend('desired','fb','ff+fb')
subplot(212)
plot(t,zd,t,x1(:,3),t,x3(:,3))
xlabel('t')
ylabel('z')
legend('desired','fb','ff+fb')

error1=yd-x3(:,1)';
error2=zd-x3(:,3)';
figure(16)
subplot(211)
plot(t,error1)
subplot(212)
plot(t,error2)
xlabel('t')
ylabel('error')
title('ff+fb trakcing error of x')
