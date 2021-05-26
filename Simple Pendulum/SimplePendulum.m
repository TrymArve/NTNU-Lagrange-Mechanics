clc;clear;close all;

GenCoords = {'th'};
S = MakeSymbolicGenCoords(GenCoords);

syms L
S.par = [L];

pos = L*arm(th);
pos_1 = pos + arm(th);
pos_2 = pos - arm(th);
S.pos = {pos,pos_1,pos_2};

S = MakeLagrange(S);
EL = S.EL;

%Fix inertia
EL = RepairInertia(EL,{'m1','m2','m3'},{'m1','I1'});
[EL,M,C,G] = ExtractMCG(EL,GenCoords);
mass = [m1;I1];
 q = [S.q{:}].';
 dq = [S.dq{:}].';
matlabFunction(M,C,G,'file','MCG','vars',{q,dq,mass,S.par});

disp('EL:')
disp(EL)

%% Reference
syms t
%thr = (sin(2*t)/(0.1*t+1))-pi/2;
%thr = 1*sin(2*t)-pi/2;
thr = pi/2+0.5*sin(2*t);
dthr = diff(thr,t);
ddthr = diff(dthr,t);
matlabFunction(thr,'file','ref','vars',t);
matlabFunction(dthr,'file','dref','vars',t);
matlabFunction(ddthr,'file','ddref','vars',t);

%% Simulate
close all;clc;
tf = 25;
L = 1;
m = 5;
I = 10;
g = 9.81;

par = [g L];
mass = [m; I];

init_state = pi/2*[1/2; 0];

Control.off = 1;
Control.PD = 1;
Control.FL = 1;

spn = Control.off+Control.PD+Control.FL;
spm = 1;
spc = 0;

color.measured = 'b-';
color.reference = 'r--';
if Control.off == 1
    % No controller:
    controller.type = "off";
    [tsim,xsim] = ode45(@(t,x) PendulumDynamics(t,x,mass,par,controller),[0 tf],init_state);
    
    %PLOT off:
    spc = spc + 1;
    subplot(spn,spm,spc)
    addleg_plot(tsim.',xsim.',color.measured,';',{'\theta','\omega';1,1});
    r = ref(tsim');
    dr = dref(tsim');
    addleg_plot(tsim',[r;dr],color.reference,';',{'\theta_{d}','\omega_{d}';1,1});
    title('No control')
end
if Control.PD == 1
    % PD controller:
    controller.PD.Kp = -L*g*m*0.5*10;
    controller.PD.Kd = -10;
    controller.type = "PD";
    [tsim,xsim] = ode45(@(t,x) PendulumDynamics(t,x,mass,par,controller),[0 tf],init_state);
    
    %PLOT PD:
    spc = spc + 1;
    subplot(spn,spm,spc)
    addleg_plot(tsim.',xsim.',color.measured,';',{'\theta','\omega';1,1});
    r = ref(tsim');
    dr = dref(tsim');
    addleg_plot(tsim',[r;dr],color.reference,';',{'\theta_{d}','\omega_{d}';1,1});
    title('PD control')
end
if Control.FL == 1
    % FL controller:
    controller.FL.Kp = 10;
    controller.FL.Kd = 1;
    controller.type = "FL";
    [tsim,xsim] = ode45(@(t,x) PendulumDynamics(t,x,mass,par,controller),[0 tf],init_state);

    %PLOT FL:
    spc = spc + 1;
    subplot(spn,spm,spc)
    addleg_plot(tsim.',xsim.',color.measured,';',{'\theta','\omega';1,1});
    r = ref(tsim');
    dr = dref(tsim');
    addleg_plot(tsim',[r;dr],color.reference,';',{'\theta_{d}','\omega_{d}';1,1});
    title('Feedback Linearization + PD control')
end


%% ANIMATE
close all;
scale = 1;

obj.ref.type = 'line';
obj.ref.a = @(q) [0;0];
obj.ref.b = @(q) (L*1.1)*arm(scale*q(2));
obj.ref.color = 'r';

obj.pend.type = {'line','ball'};
obj.pend.a = @(q) [0;0];
obj.pend.b = @(q) L*arm(scale*q(1));
obj.pend.color = 'c';
obj.pend.c = obj.pend.b;
obj.pend.r = .1;

lift = 0.5;
shift = 0;
config.axis = L*[-1.25 1.25 -1 1]*1.1 + [shift shift lift lift];
config.simspeed = 2;
config.tf = 20;
config.grid = 'on';
Animate(tsim,[xsim(:,1) r'],obj,config)

%% Functions

function[dstate] = PendulumDynamics(t, state, mass, par, controller)

q  = state(1);
dq = state(2);

q_ref = ref(t);
dq_ref = dref(t);
ddq_ref = ddref(t);

g = par(1);
L = par(2);
m = mass(1);
I = mass(2);

%CONTROL:
c = controller.type;

if c == "off"
    u = 0;
elseif c == "PD"
    Kp = controller.PD.Kp;
    Kd = controller.PD.Kd;
    ep  = q  - q_ref;
    ed  = dq - dq_ref; 
    u = Kp*ep + Kd*ed;
elseif c == "FL"
    Kp = controller.FL.Kp;
    Kd = controller.FL.Kd;
    ep  = q  - q_ref;
    ed  = dq - dq_ref; 
    v = ddq_ref - Kp*ep - Kd*ed; %Assumes q_ref is C^2(continuousl differentiable twice)
    Lin = L*g*m*cos(q); % Linearizing term
    u = (m*L^2 + I)*v + Lin;
end

U = [u];

[M,C,G] = MCG(q,dq,mass,par);
ddq = M\(U - C*dq - G);


dstate = [dq;ddq];

end