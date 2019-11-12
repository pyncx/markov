function kcnq1iv

% I-V Curve

vact=[-70:10:40]';  % Potentials to depolarize to

erev=-90;           % Reversal Potential (mV)
y0=[1 zeros(1,20)]; % Initial Conditions

[y1 y2]=ode23s(@kcnq1rates,[0 1e7],y0,[],-80);  % Steady-State at holding potential

yp=y2(length(y2),:);

for i=1:length(vact)
    [y1 y2]=ode23s(@kcnq1rates,[0:2000],yp,[],vact(i));    % Activation 
    [y3 y4]=ode23s(@kcnq1rates,[0:2000],y2(length(y2),:),[],-70);  %  Deactivation
    Iact(:,i)=sum(y2(:,16:20),2)*(vact(i)-erev);                % Find Current  
    Ideact(:,i)=sum(y4(:,16:20),2)*(-70-erev);
end

ttp=[y1;y3+max(y1)];                    % Concatenate time vectors
Itp=[Iact;Ideact]/max(max(Iact));       % Concatenate currents

figure                                  % Plot
plot(ttp,Itp)
title('KCNQ1 IVcurve')
xlabel('Time (ms)')
ylabel('Current (Arbitrary Units)')

function y=kcnq1rates(t,x,v)

F=96485;
T=298;
R=8314;

alpha=9.57e-004*exp(v*F/(R*T)*1.98e-001);
beta=5.00e-005*exp(v*F/(R*T)*-3.33e-002);
gamma=3.77e-002*exp(v*F/(R*T)*3.33e-002);
delta=4.77e-002*exp(v*F/(R*T)*-4.77e-001);
eta=7.98e-002;
zeta=2.06e-002*exp(v*F/(R*T)*-1.04e+000);
epsilon=2.54e-002*exp(v*F/(R*T)*6.46e-002);
theta=1.78e-002*exp(v*F/(R*T)*-5.32e-001);
lambda=2.32e-002*exp(v*F/(R*T)*1.21e-001);
mu=6.19e-002*exp(v*F/(R*T)*-9.46e-002);
xe=8.74e-001;
xz=3.29e-001;

c1=x(1);
c2=x(2);
c3=x(3);
c4=x(4);
c5=x(5);
c6=x(6);
c7=x(7);
c8=x(8);
c9=x(9);
c10=x(10);
c11=x(11);
c12=x(12);
c13=x(13);
c14=x(14);
a=x(15);
o1=x(16);
o2=x(17);
o3=x(18);
o4=x(19);
o5=x(20);
i1=x(21);

di1dt = lambda*o5-mu*i1;
do5dt = epsilon*xe^3*o4+mu*i1-(4*theta*xz^3+lambda)*o5;
do4dt = 2*epsilon*xe^2*o3+4*theta*xz^3*o5-(3*theta*xz^2+epsilon*xe^3)*o4;
do3dt = 3*epsilon*xe*o2+3*theta*xz^2*o4-(2*theta*xz+2*epsilon*xe^2)*o3;
do2dt = 4*epsilon*o1+2*theta*xz*o3-(theta+3*epsilon*xe)*o2;
do1dt = eta*a+theta*o2-(4*epsilon+zeta)*o1;
dadt = gamma*c14+zeta*o1-(4*delta+eta)*a;
dc14dt = alpha*c13+4*delta*a+2*gamma*c12-(beta+3*delta+gamma)*c14;
dc13dt = beta*c14+gamma*c11-(alpha+3*delta)*c13;
dc12dt = alpha*c11+3*delta*c14+3*gamma*c9-(2*beta+2*delta+2*gamma)*c12;
dc11dt = 2*alpha*c10+2*beta*c12+2*gamma*c8+3*delta*c13-(beta+alpha+gamma+2*delta)*c11;
dc10dt = beta*c11+gamma*c7-(2*alpha+2*delta)*c10;
dc9dt = alpha*c8+2*delta*c12+4*gamma*c5-(3*beta+delta+3*gamma)*c9;
dc8dt = 2*alpha*c7+3*beta*c9+3*gamma*c4+2*delta*c11-(2*beta+alpha+2*gamma+delta)*c8;
dc7dt = 3*alpha*c6+2*beta*c8+2*gamma*c3+2*delta*c10-(beta+2*alpha+delta+gamma)*c7;
dc6dt = beta*c7+gamma*c2-(3*alpha+delta)*c6;
dc5dt = alpha*c4+delta*c9-(4*beta+4*gamma)*c5;
dc4dt = 2*alpha*c3+4*beta*c5+delta*c8-(3*beta+alpha+3*gamma)*c4;
dc3dt = 3*alpha*c2+3*beta*c4+delta*c7-(2*beta+2*alpha+2*gamma)*c3;
dc2dt = 4*alpha*c1+2*beta*c3+delta*c6-(beta+3*alpha+gamma)*c2;
dc1dt = beta*c2-4*alpha*c1;

y(1,1)=dc1dt;
y(2,1)=dc2dt;
y(3,1)=dc3dt;
y(4,1)=dc4dt;
y(5,1)=dc5dt;
y(6,1)=dc6dt;
y(7,1)=dc7dt;
y(8,1)=dc8dt;
y(9,1)=dc9dt;
y(10,1)=dc10dt;
y(11,1)=dc11dt;
y(12,1)=dc12dt;
y(13,1)=dc13dt;
y(14,1)=dc14dt;
y(15,1)=dadt;
y(16,1)=do1dt;
y(17,1)=do2dt;
y(18,1)=do3dt;
y(19,1)=do4dt;
y(20,1)=do5dt;
y(21,1)=di1dt;

if abs(sum(y))>1e-6
    error('sum of derivatives greater than 0')
end
