clc;clear;close all;

GenCoords = {'q1','q2'};
Parameters = {'L_cm','L'};
S = MakeSymbolicGenCoords(GenCoords,Parameters);
q1 = q1 + pi/2; % q1 is zero when pointing upwards
pos_pend = L_cm*arm(q1);
pos_pend_1 = pos_pend + arm(q1);
pos_pend_2 = pos_pend - arm(q1);
pos_wheel = L*arm(q1);
pos_wheel_1 = pos_wheel + arm(q2);
pos_wheel_2 = pos_wheel - arm(q2);
S.pos = {pos_pend,pos_pend_1,pos_pend_2,pos_wheel,pos_wheel_1,pos_wheel_2};

% Define your own mass-names(not necessary, if "MakeLagrange" doesn't find an "S.mass", it will simply make it own masses)
syms mp mp1 mp2 mw mw1 mw2
S.mass = {mp, mp1, mp2, mw, mw1, mw2};

S.DisableExtension = 1;  % We don't need to make the MCG-function before we have conveted it to include the inertia
S = MakeLagrange(S);

%% Convert EL-equation to include the inertia, rather than the two stretched out masses:
syms Iw Ip
S.EL = RepairInertia(S.EL,{'mw','mw1','mw2'},{'mw','Iw'}); %Repair inertia of the pole
S.EL = RepairInertia(S.EL,{'mp','mp1','mp2'},{'mp','Ip'}); %Repair the interia of the wheel

S.mass = {mw,Iw,mp,Ip}; % Remember to define new masses/inertias in your system struct
S.mass_ = [S.mass{:}].';

disp('EL(fixed inertia):')
disp(S.EL)

%% Extract MCG-form (manually, this time)

%Normally, MakeLagrange does this for us, but, since we changed the
%EL-equation afterwards, we need to do this with the new EL-expression

[EL,M,C,G] = ExtractMCG(S.EL,GenCoords); 
matlabFunction(M,C,G,'file','MCG','vars',{S.q_,S.dq_,S.mass_,S.par});

%% Define masses and parameters
m1 = 1;
I1 = 0.5;
m2 = 2;
I2 = 3;
mass = [m1 I1 m2 I2]';

g = 9.81;
L_cm = 2;
L = 4;
parameters = [g L_cm L]';

%% Configure Simulation

init_coords = [0.3 ; -0.2];
init_vel = [0;0];
init_state = [init_coords;init_vel];

%% controller 1: Actuation on pole and wheel, pole-placement-control
clear("controller")
tf = 5;
S.B = eye(2);
controller.type = "KF";
S = LinearizeEL(S,[0;0],[0;0],mass,parameters); % Linearize about upwards equilibrium
poles = -[1 2 3 4];
controller.KF.K = place(double(S.LA),double(S.LB),poles);

%% Controller 2: Actuation on wheel, PD control (no measurement on wheel angle)
clear("controller")
tf = 10;
controller.type = "PD";
controller.PD.C = [1 0]; %measure pole angle
controller.PD.Kp = 1000;
controller.PD.Kd = 100;

%% Controller 3: Actuation on wheel, PD control, stabilizing BOTH coordinates
clear("controller")
tf = 50;
controller.type = "PD";
controller.PD.C = eye(2); %measure both angles
controller.PD.Kp = [1000 -1];
controller.PD.Kd = [100 -1];
values.umax = 100;
%% Simulate

Simulate_EL

%% Animate
close all;
% Define object to animate:
clear("obj")

X = xsim(:,1:end/2);
X(:,1) = X(:,1) + pi/2;

%pole
obj.pole.type = 'line';
obj.pole.a = @(q) [0;0];
obj.pole.b = @(q) L*arm(q(1));
obj.pole.color = 'b';

%wheel
radius_w = 0.5;
obj.wheel.type = {'line','ball'};
obj.wheel.c = @(q) L*arm(q(1));
obj.wheel.r = radius_w;
obj.wheel.a = obj.wheel.c;
obj.wheel.b = @(q) L*arm(q(1)) + radius_w*arm(q(1)+q(2));

if controller.type == "KF"
%Torque on pole
obj.torque1.type = 'moment';
obj.torque1.target = obj.pole.a;
obj.torque1.magnitude = @(q) q(values.nq+1)*0.1;
obj.torque1.color = GetColorCode('r');
obj.torque1.color2 = GetColorCode('b',1.1);
end

%Torque on wheel
obj.torque2.type = 'moment';
obj.torque2.target = obj.wheel.c;
obj.torque2.magnitude = @(q) q(end)*0.1;
obj.torque2.color = GetColorCode('r');
obj.torque2.color2 = GetColorCode('b',1.1);



% Configure animation settings:
config.simspeed = 2;
config.axis = [-4.6 4.6 -3 5]*1.5;
config.tf = tf;
config.enterToStart = 0;
Animate(tsim,[X usim],obj,config);
