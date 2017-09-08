clc 
clear all
close all
format short e
%step 1  state space model
%system parameters
lamda=1;eplsion=0.1;M=1;J=1;g=9.81;
%mass matrix
A=[0 1 0 0 0 0;
    0 0 0 0 -g 0;
    0 0 0 1 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 1;
    0 0 0 0 0 0];
B=[0 0;0 eplsion/M;0 0;1/M 0;0 0;0 lamda/J];
C=[1 0 0 0 0 0;0 0 1 0 0 0];
C1=C(1,:);C2=C(2,:);
D=0;

%step 2----the desired output trajectory
T=1;%T is the time period in seconds
om=2*pi/T;
delt=0.001;
Ymax=1;
Ay=Ymax*2*pi/(T^2);

t=0:delt:T;
yd=(Ay/om)*t-(Ay/om/om)*sin(om*t);
ydd=(Ay/om)*(1-cos(om*t));
yd2d=Ay*sin(om*t);

figure(2)
subplot(311),plot(t,yd)
ylabel('position y_d')
subplot(312),plot(t,ydd)
ylabel('velocity')
subplot(313),plot(t,yd2d)
ylabel('acceleration')
xlabel('time')

%add zeros/constants at the beginning and final
Npre=4;Npost=Npre;
tpre=0:delt:(Npre*T-delt);
tpost=0:delt:(Npost*T-delt);

t=0:delt:(1+Npre+Npost)*T;

yd=[zeros(size(tpre)) yd Ymax*ones(size(tpost))];
ydd=[zeros(size(tpre)) ydd zeros(size(tpost))];
yd2d=[zeros(size(tpre)) yd2d zeros(size(tpost))];

figure(1)
subplot(311),plot(t,yd)
ylabel('position y_d')
subplot(312),plot(t,ydd)
ylabel('velocity')
subplot(313),plot(t,yd2d)
ylabel('acceleration')
xlabel('time')

%step3---internal dynamics
Tt=[C1;C1*A;C2;C2*A];
Tb=[0 0 0 0 1 0;0 0 0 0 0 1];
T=[Tt;Tb];
invT=inv(T);
invTleft=invT(:,1:4);
invTright=invT(:,4:6);
Ay=[C1*A*A;C2*A*A];
By=[C1*A*B;C1*A*B];

Ainv=Tb*A*invTright-Tb*B*inv(By)*Ay*invTright;
Binv=[(Tb*A*invTleft-Tb*B*inv(By)*Ay*invTleft) Tb*B*inv(By)];

[V,poles_internal_dynamics]=eig(Ainv);
poles_internal_dynamics

Tsplit=V;
invTsplit=inv(Tsplit);

A_split=invTsplit*Ainv*Tsplit;
B_split=invTsplit*Binv;
A_s=A_split(1,1);B_s=B_split(1,:);
A_u=A_split(2,2);B_u=B_split(2,:);

%step 4 solve the internal dynamics to find inverse 
capYd=[yd;ydd;yd2d];
figure(50);plot(t,capYd)
legend('y','dy/dt','d^2y/dt^2')
xlabel('time')
ylabel('Y_d')

sys_s=ss(A_s,B_s,[1],[0 0 0]);

sys_u_backward=ss(-A_u,-B_u,[1],[0 0 0]);

%define initial and final states as equilbrium states
X_initial=[0 0 0 0 0 0]';
X_final=[Ymax 0 Ymax 0 0 0]';
%initial conditions in transformed co-ordinates
Zeta_eta_initial=T*X_initial;
Zeta_eta_final=T*X_final;
%initial conditions of internal states
eta_initial=Zeta_eta_initial(3:4,1);
eta_final=Zeta_eta_final(3:4,1);
%initial conditions in split co-ordinates
eta_split_initial=invTsplit*eta_initial;
eta_split_final=invTsplit*eta_final;
eta_s_initial=eta_split_initial(1,1);
eta_u_initial=eta_split_initial(2,1);
eta_s_final=eta_split_final(1,1);
eta_u_final=eta_split_final(2,1);

%solve the stable internal dynamics
[eta_s,time,x_s]=lsim(sys_s,capYd,t,eta_s_initial);

figure(51);plot(t,eta_s)
xlabel('t')
ylabel('stable part \eta_{s,ref}')
grid

%use this to simulate error in initial condition
[eta_s2,time,x_s]=lsim(sys_s,capYd,t,-0.2);
figure(71);plot(t,eta_s,'b:',t,eta_s2,'r')
legend('without IC error','with IC error')
xlabel('t')
ylabel('stable part \eta_{s,ref}')
% axis([0 max(t) -0.8 0.1])]
%eta_s=eta_s2; %uncomment this if you want to simulate

%solve the unstable internal dynamics
Capydbackward=fliplr(capYd);
figure(52);plot(t,Capydbackward)
legend('y','dy/dt','d^2y/dt^2')
xlabel('t');ylabel('Y_{d,backward}')

[eta_u_backward,time,x_u]=lsim(sys_u_backward,Capydbackward,t,eta_u_final);

figure(53);plot(t,eta_u_backward)
xlabel('t');ylabel('\eta_u backward')
grid

%flip the unstable solution to the correct time directly
eta_u=fliplr(eta_u_backward');

figure(54);plot(t,eta_u)
xlabel('t');ylabel('\eta_{u,ref}')
grid

%use this to simulate error in final condition
[eta_u2,time,x_u]=lsim(sys_u_backward,Capydbackward,t,0);
eta_u2=fliplr(eta_u2');
figure(72);plot(t,eta_u,'b:',t,eta_u2,'r')
legend('without IC error','with IC error')
xlabel('t')
ylabel('unstable part \eta_{u,ref}')

%try to simulate forward in time
sys_uf=ss(A_u,B_u,[1],[0 0 0]);
[eta_uf,time,x_u]=lsim(sys_uf,capYd,t,0);
%[eta_uf,time,x_u]=lsim(sys_uf,capYd,t,-1.3559e-09);
figure(73);plot(t,eta_uf)
xlabel('t')
ylabel('\eta_{u,forward}')
grid

%plot the results for the split internal dynamics
figure(3)
subplot(211),plot(t,eta_s)
ylabel('\eta_s')
subplot(212),plot(t,eta_u)
ylabel('\eta_u')
xlabel('t')

%find the origianl internal dynamics eta_ref
eta_ref=Tsplit*[eta_s';eta_u];
figure(4);
plot(t,eta_ref)
ylabel('\eta_{ref}')

%find the reference state trajectory xref
zeta_ref=[yd;ydd];
xref=invT*[zeta_ref;eta_ref];

%the inverse input
Uff=(yd2d-C*A*A*xref)/(C*A*B);

figure(55)
plot(t,Uff)
xlabel('t')
ylabel('Inverse input u_{ff}')
grid


%design feedback controller
%use LQR method to find Kfb
%an alternative is to use "place"
%try changing R and Q
R=0.01;Q=C'*C/.10;[Kfb,S,E]=lqr(sys,Q,R);

Acl=A-B*Kfb;
Ut=Uff+Kfb*xref;
sysCL=ss(Acl,B,C,D);

x0=[0 0 0 0 0 0]';

[ycl,time,xcl]=lsim(sysCL,Ut,t,x0);

%the actual input to the system
Ufb=-Kfb*(xcl'-xref);
Utotal=Uff+Ufb;

%plot the outputs and the inputs
figure(5)
plot(t,ycl,'r',t,yd,'g')
legend('y_{cl}','y_d')
xlabel('t')
ylabel('output')

figure(6)
plot(t,Uff,'r',t,Ufb,'g',t,Utotal,'b')
legend('u_{ff}','u_{fb}','u')
xlabel('t')
ylabel('inputs')
%check the size of feedback input
%compare to the feedforward input

U_fbff=Utotal;
y_fbff=ycl;

%simulate the response to feedback alone
r=(10/0.3154)*yd;%scale r to get zero steady state

[ycl,time,xcl]=lsim(sysCL,r,t);
%the actual input to the system
Ufb=-Kfb*(xcl');
Utotal=r+Ufb;
U_fb=Utotal;
y_fb=ycl;

%comparison of tracking error, input and output in both direction
figure(7)
subplot(311),plot(t,U_fbff,'g',t,U_fb,'r:')
legend('u_{fb+ff}','u_{fb}')
xlabel('t')
ylabel('inputs')

subplot(312),plot(t,yd,'b',t,y_fbff,'g',t,y_fb,'r:')
legend('yd','y_{fb+ff}','y_{fb}')
xlabel('t')
ylabel('outputs')

subplot(313),plot(t,yd-y_fbff','g',t,yd-y_fb','r:')
legend('e_{fb+ff}','e_{fb}')
xlabel('t')
ylabel('eror y-y_d')



