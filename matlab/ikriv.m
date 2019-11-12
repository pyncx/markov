function ikriv

% I-V Curve

k_o=4;      % mM
erev=-94;   % Reversal potential (mV)
x0=[1 0 0 0 0];

ivv=[-30:10:40]';   % Potentials to depolarize to

[y1 y2]=ode23s(@ikrrates,[0 10000],x0,[],-40,k_o);    % Steady-state at holding potential

yp=y2(length(y2),:);

for i=1:length(ivv)
    
    [y1 y2]=ode23s(@ikrrates,[0:550],yp,[],ivv(i),k_o);
    Iact(:,i)=y2(:,4)*(ivv(i)-erev);
    [y3 y4]=ode23s(@ikrrates,[0:500],y2(length(y2),:),[],-40,k_o);
    Itail(:,i)=y4(:,4)*(-40-erev);
    
end

figure  %   plot

title('IKr IVcurve')
xlabel('Time (ms)')
ylabel('Current (Arbitrary Units)')
plot([y1;y3+max(y1)],[Iact;Itail])

function y=ikrrates(t,x,v,k_o)

F=96485;
T=310;
R=8314;

alpha2=1.31e-2*exp(1.48*v*F/(R*T));
alpha1=2.17;
alpha=3.02e-2*exp(1.48*v*F/(R*T));
beta=2.90e-3*exp(-9.78e-1*v*F/(R*T));
beta1=1.08;
beta2=3.30e-3*exp(-5.77e-1*v*F/(R*T));
alphai=5.45e-1*exp(-8.17e-1*v*F/(R*T))*4.5/k_o;
betai=8.20e-1*exp(5.04e-1*v*F/(R*T))*(4.5/k_o)^0.3;
mu=(alphai*beta2)/betai;

c3=x(1);
c2=x(2);
c1=x(3);
o=x(4);
i=x(5);

dc3dt=beta*c2-alpha*c3;
dc2dt=beta1*c1+alpha*c3-(alpha1+beta)*c2;
dc1dt=alpha1*c2+beta2*o+mu*i-(beta1+2*alpha2)*c1;
dodt=alphai*i+alpha2*c1-(betai+beta2)*o;
didt=alpha2*c1+betai*o-(mu+alphai)*i;

y(1,1)=dc3dt;
y(2,1)=dc2dt;
y(3,1)=dc1dt;
y(4,1)=dodt;
y(5,1)=didt;

if abs(sum(y))>1e-6
    error('sum of derivatives greater than 0')
end
