

%% Define Dynamics
clc;clear;close all;

S = MakeSymbolicGenCoords({'th'},{'L'});

S.pos = {arm(th+pi/2)*L};  % Zero is upwards

S = MakeLagrange(S);


%% Define parameter values
g = 9.81;
L = 1;
m1 = 2;

parameters = [g;L];
mass = [m1];
%% Configure Simulation

tf = 5;
init_state = [2;0]; %[th;dth]

%% Configure Controller

clear("controller")
controller.type = "PD";
%controller.type = "off";  %Turns off controller
controller.PD.Kp = 200;
controller.PD.Kd = 20;
controller.PD.C = [1]; % send the error through a measurement matrix before controlling (u = Kp*C*(ref - q) + Kd*c*(dref - dq))

values.umax = 30; %Maximum torque that the motor can apply

%% Simulate
Simulate_EL

%% Animate

close all;
scale = 1;

obj.ref.type = {'line','ball'};
obj.ref.a = @(q) [0;0];
obj.ref.b = @(q) L*arm(q(3)+pi/2);
obj.ref.color = GetColorCode('i',1);
obj.ref.c = obj.ref.b;
obj.ref.r = .08;

obj.pend.type = {'line','ball'};
obj.pend.a = @(q) [0;0];
obj.pend.b = @(q) L*arm(q(1)+pi/2);
obj.pend.color = GetColorCode('g',1);
obj.pend.c = obj.pend.b;
obj.pend.r = .1;

obj.input.type = 'arrow';
obj.input.tail = @(q) obj.pend.b(q);
obj.input.head = @(q) obj.pend.b(q) + 0.1*q(2)*arm(q(1)+pi);
obj.input.color = GetColorCode('r');


% CONFIGURE ANIMATION:
formatRatio = 5/4;
formatRatio = 5/4*1.55;
formatRatio = 5/4*0.75;
lift = 0;
shift = 0;
height = 3*L;
width = height*formatRatio;

config.axis = [-width/2 width/2 -height/2 height/2] + [shift shift lift lift];
config.simspeed = 1;
config.tf = tf;
config.grid = 'on';
config.enterToStart = 1;
Animate(tsim,[xsim(:,1) usim Ref],obj,config);










