
%% Define Dynamics
clc;clear;close all;

S = MakeSymbolicGenCoords({'th','x'},{'L'});

pos_cart     = [x;0];
pos_pendulum = pos_cart + arm(th+pi/2)*L;

S.pos = {pos_cart,pos_pendulum};

S = MakeLagrange(S);

%% Give parameter values
g = 9.81;
L = 1;
m1 = 1;
m2 = 1;


parameters = [g;L];
mass = [m1;m2];

%% Configure Simulation

tf = 10; % Duration

initial_coordinates = [1; -0.5];  %[th ,x]
initial_velocities  = [0;0.2];    %[dth,dx]
init_state = [initial_coordinates; initial_velocities];

%% Simulate
Simulate_EL

%% Animate

close all;clear("obj");clear("config");

box_width = 0.5;
box_angle = pi/6;

%Cart:
obj.cart.type = 'box';
obj.cart.def = '2C';
obj.cart.B1 = @(q)[q(2);0] + box_width*arm(pi+box_angle);
obj.cart.B2 = @(q)[q(2);0] + box_width*arm(box_angle);
obj.cart.color = GetColorCode('b',1.2);

%Pendulum:
obj.pend.type = {'line','ball'};
obj.pend.a = @(q) [q(2);0];
obj.pend.b = @(q) [q(2);0] + L*arm(q(1)+pi/2);
obj.pend.color = GetColorCode('g');
obj.pend.c = obj.pend.b;
obj.pend.r = .1;


% CONFIGURE ANIMATION:

config.position = [1 1 1920 1080];
config.framecenter = [0 -0.2];
config.frameheight = 3*L; 
config.figureheight = 1080;
config.simspeed = 1;
config.tf = tf;
config.grid = 'on';
config.enterToStart = 0;

% Run animation
Animate(tsim,xsim(:,1:end/2),obj,config);



























