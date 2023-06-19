function RCTA_PID
clc
tint=[0,100];
x0=[2.7714,335.2770,338.2331,0,0,0,0,0];
Tref=330;
[t,x]=ode23tb(@(t,x)Sistema(t,x,Tref,0),tint,x0);
n=length(t);
vj=zeros(n,1);
for i=1:1:n
    vj(i,1)=Sistema(t(i,1),x(i,:),Tref,1);
end

A={'$C_A \qquad \left(\frac{mol}{m^3}\right) $',...
    '$T \qquad \left(K\right) $',...
    '$T_j \qquad \left(K\right) $'};
B= {'$t\qquad \left(s\right) $'};
C= {'$v_J\qquad \left(\frac{L}{s}\right) $'};
for i=1:1:3
    figure(i)
    subplot(2,1,1)
    plot(t,x(:,i),'m','LineWidth',2)
    grid on, grid minor
    xlabel(B,'FontSize',14,'Interpreter',...
        'latex')
    ylabel(A{1,i},'FontSize',14,'Interpreter',...
        'latex')
    hold on
    subplot(2,1,2)
    plot(t,vj,'m','LineWidth',2)
    grid on, grid minor
    xlabel(B,'FontSize',14,'Interpreter',...
        'latex')
    ylabel(C,'FontSize',14,'Interpreter',...
        'latex')
    hold on
end


end
function dx=Sistema(~,x,Tref,w)
%% Variables de estado
CA=x(1); T=x(2); Tj=x(3);
%%Variables del sensor
H=x(4);
%% Variables de la valvula
Va=x(5); va=x(6);
%% Variables del controlador
Z1=x(7); Z2=x(8);
%% Parámetros
v=0.1;    %m3/s
vjs=0.25;     %m3/s
CA0=40;    %mol/m3
T0=323.15; %K
Tj0=353.15; %K
V=1; %m3
rho=1000; % kg/m3
Cp=4184; %J/kgK
Vj=1;     %m3
Cpj=1514; %J/kgK
rhoj=658; %kg/m3
UA=1256800; %J/K s
DHrxn=-365000;%J/mol
k1=1; T1=50+273.15;
k2=2; T2=80+273.15;
%% Ecuaciones auxiliares
ER=-log(k2/k1)/(1/T2-1/T1);
k0=k1*exp(ER/T1);
k=k0*exp(-ER/T);
%% Parametros del sensor
Km=(20-4)/(750-0);
taum=0.05;
%% Parametros de la valvula
Kv=(0.5-0)/(20-4);
deltav=2;
tauv=0.3;
%% Parametros del controlador
Kc=402.373;
tauI=3.844;
tauD=0.961;
a=0.15/tauD;
%% Salida
Ts=335.2770;
y=T-Ts;
%% Error
R=Tref-Ts;
E=Km*R-H;
%% Controlador
C=Kc*((1+1/a)*E+Z1+Z2);
%% Entrada
vj=vjs+Va;
%% Limitaciones fisicas
if vj<0
    vj=0;
end
if vj>0.5
    vj=0.5;
end
if w==0
%% Ecuaciones diferenciales
dCA=v*(CA0-CA)/V-k*CA;
dT=v*(T0-T)/V-UA*(T-Tj)/(rho*V*Cp)...
    -DHrxn*k*CA/(rho*Cp);
dTj=vj*(Tj0-Tj)/Vj+UA*(T-Tj)/(rhoj*Vj*Cpj);
dH=(Km*y-H)/taum;
dVa=va;
dva=-Va/tauv^2-2*deltav*va/tauv+Kv*C/tauv^2;
dZ1=E/tauI;
dZ2=-(E/a+Z2)/0.15;
dx=[dCA;dT;dTj;dH;dVa;dva;dZ1;dZ2];
else
    dx=vj*1000;
end

end