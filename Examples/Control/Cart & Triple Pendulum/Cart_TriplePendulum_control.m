
%% Define System (2D)
clc;clear;close all;

S = MakeSymbolicGenCoords({'x','th','phi','xi'},{'L1','L2','L3'}); 

pos_cart   = [x;0];
pos_pend_1 = pos_cart   + L1*arm(th  + pi/2);
pos_pend_2 = pos_pend_1 + L2*arm(phi + pi/2);
pos_pend_3 = pos_pend_2 + L3*arm(xi  + pi/2);
S.pos = {pos_cart,pos_pend_1,pos_pend_2,pos_pend_3};

S = MakeLagrange(S);

%% Give values to the various system parameters and masses:

% Parameters:
g = 9.81;
L1 = 0.5;
L2 = 1;
L3 = 2;
parameters = [g;L1;L2;L3];

% Masses:
m1 = 1;
m2 = 1.5;
m3 = 2;
m4 = 0.5;
mass = [m1;m2;m3;m4];

%% Configure the Simulation:

tf = 15;

% Initial State for the simulation:
initial_coordinates = [0 0 0 0].' + pi;  %[x; th; phi; xi]
initial_velocities  = [0 0 0 0].';  %[dx; dth; dphi; dxi]
init_state = [initial_coordinates; initial_velocities];

% Maximum input value(saturation level)
values.umax = 75;

%% Configure Controller:

% Only input the cart:
controller.B = [1 0 0 0].';

% Select what controller to use:
 controller.type = "PD";

% PD:
controller.PD.C  = [1 0 0 0]; % Measure output (y = C*q) --> (e = C*(q - ref))
controller.PD.Kp = 100; %Proportional gain
controller.PD.Kd = 10; %Derivative gain

%% Set the Reference Trajectory:
freq = 1;
amp = 1;
controller.ref   = @(t) [amp*sin(freq.*t);         0.*t; 0.*t; 0.*t]; % Add ".*t" so that all elements are of the same size when inputing a timeseries
controller.dref  = @(t) [amp*freq*cos(freq.*t);    0.*t; 0.*t; 0.*t]; 
controller.ddref = @(t) [-amp*freq^2*sin(freq.*t); 0.*t; 0.*t; 0.*t]; 


%% SIMULATION  
Simulate_EL

%% Animate

close all;
scale = 1;
clear("obj")

AnimationValues = [xsim(:,1:end/2) usim Ref]; % All coordinates that should be send to the "Animate" function. i.e. all coordinates that are needed to plot the various object in the animation.
startU   = size(xsim(1,1:end/2),2);
startRef = startU + size(usim,2);

CartHeight = 0.25;
CartDiagonal = 0.5;

DisableReference = 1;

%Ref Cart:
obj.ref_cart.type = 'box';
obj.ref_cart.def = "CD";
obj.ref_cart.B1 = @(q) [q(startRef+1); 0];
obj.ref_cart.B2 = @(q) CartDiagonal;
obj.ref_cart.color = GetColorCode('i',1.1);
obj.ref_cart.disable = DisableReference*0;
%Ref Pendulum 1:
obj.ref_1.type = {'line','ball'};
obj.ref_1.a = @(q) obj.ref_cart.B1(q)+[0;CartHeight/2];
obj.ref_1.b = @(q) obj.ref_1.a(q) + L1*arm(q(startRef+2)+pi/2);
obj.ref_1.color = GetColorCode('i',1.1);
obj.ref_1.c = obj.ref_1.b;
obj.ref_1.r = .05;
obj.ref_1.disable = DisableReference;
%Ref Pendulum 2:
obj.ref_2.type = {'line','ball'};
obj.ref_2.a = @(q) obj.ref_1.b(q);
obj.ref_2.b = @(q) obj.ref_2.a(q) + L2*arm(q(startRef+3)+pi/2);
obj.ref_2.color = GetColorCode('i',1.1);
obj.ref_2.c = obj.ref_2.b;
obj.ref_2.r = .05;
obj.ref_2.disable = DisableReference;
%Ref Pendulum 3:
obj.ref_3.type = {'line','ball'};
obj.ref_3.a = @(q) obj.ref_2.b(q);
obj.ref_3.b = @(q) obj.ref_3.a(q) + L2*arm(q(startRef+4)+pi/2);
obj.ref_3.color = GetColorCode('i',1.1);
obj.ref_3.c = obj.ref_3.b;
obj.ref_3.r = .05;
obj.ref_3.disable = DisableReference;

%Cart:
obj.cart.type = 'box';
obj.cart.def = "CD";
obj.cart.B1 = @(q) [q(1); 0];
obj.cart.B2 = @(q) 0.5;
obj.cart.color = GetColorCode('o',0.8);

%Pendulum 1:
obj.pend_1.type = {'line','ball'};
obj.pend_1.a = @(q) obj.cart.B1(q) + [0;CartHeight/2];
obj.pend_1.b = @(q) obj.pend_1.a(q) + L1*arm(q(2)+pi/2);
obj.pend_1.color = GetColorCode('g',0.8);
obj.pend_1.c = obj.pend_1.b;
obj.pend_1.r = .05;

%Pendulum 2:
obj.pend_2.type = {'line','ball'};
obj.pend_2.a = @(q) obj.pend_1.b(q);
obj.pend_2.b = @(q) obj.pend_2.a(q) + L2*arm(q(3)+pi/2);
obj.pend_2.color = GetColorCode('b',0.8);
obj.pend_2.c = obj.pend_2.b;
obj.pend_2.r = .05;

%Pendulum 3:
obj.pend_3.type = {'line','ball'};
obj.pend_3.a = @(q) obj.pend_2.b(q);
obj.pend_3.b = @(q) obj.pend_3.a(q) + L2*arm(q(4)+pi/2);
obj.pend_3.color = GetColorCode('p',0.8);
obj.pend_3.c = obj.pend_3.b;
obj.pend_3.r = .05;

%Input torque on Cart:
obj.input_cart.type = 'arrow';
obj.input_cart.head = @(q) obj.cart.B1(q) + (1 - 2*(q(startU + 1) > 0))*[1;0].*cos(pi/6)*CartDiagonal/2;
obj.input_cart.tail = @(q) obj.input_cart.head(q) + [-1;0]*q(startU + 1)/10;
obj.input_cart.color = GetColorCode('r',0.9);

% %Input torque on Pundulum 1:
% obj.input_pend_1.type = 'arrow';
% obj.input_pend_1.tail = @(q) obj.pend_1.b(q);
% obj.input_pend_1.head = @(q) obj.pend_1.b(q) + 0.1*q(startU+2)*arm(q(2)+pi);
% obj.input_pend_1.color = GetColorCode('r',0.9);
% 
% %Input torque on Pundulum 2:
% obj.input_pend_2.type = 'arrow';
% obj.input_pend_2.tail = @(q) obj.pend_2.b(q);
% obj.input_pend_2.head = @(q) obj.pend_2.b(q) + 0.1*q(startU+3)*arm(q(3)+pi);
% obj.input_pend_2.color = GetColorCode('r',0.9);
% 
% %Input torque on Pundulum 2:
% obj.input_pend_3.type = 'arrow';
% obj.input_pend_3.tail = @(q) obj.pend_3.b(q);
% obj.input_pend_3.head = @(q) obj.pend_3.b(q) + 0.1*q(startU+4)*arm(q(4)+pi);
% obj.input_pend_3.color = GetColorCode('r',0.9);


% CONFIGURE ANIMATION:
config.frameheight = 2*(L1+L2+L3);
config.aspect = 5/4; %Change this to whatever suits your pop-up figure
config.simspeed = 1;
config.tf = tf;
config.grid = 'on';
config.enterToStart = 1;
Animate(tsim,AnimationValues,obj,config);



