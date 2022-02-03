
%% Define Dynamics
clc;clear;close all;

S = MakeSymbolicGenCoords({'th','phi','xi'},{'L1','L2','L3'});

pos_pend_1 =              L1*arm(th+pi/2);
pos_pend_2 = pos_pend_1 + L2*arm(phi+pi/2);
pos_pend_3 = pos_pend_2 + L3*arm(xi+pi/2);

S.pos = {pos_pend_1,pos_pend_2,pos_pend_3};

S = MakeLagrange(S);


%% Resulting Expressions
clc;
disp('Euler-Lagrange equations of motion:')
disp(S.EL)

disp('Kinetic energies of masses:')
disp('m1:')
disp(S.Ki{1})
disp('m2:')
disp(S.Ki{2})
disp('m3:')
disp(S.Ki{3})
disp('total kinetic energy:')
disp(S.K)

disp('Potential energies of masses:')
disp('m1:')
disp(S.Pi{1})
disp('m2:')
disp(S.Pi{2})
disp('m3:')
disp(S.Pi{3})
disp('total potential energy:')
disp(S.P)


%% Give parameter values

g = 9.81;
L1 = 1;
L2 = 1;
L3 = 1;
m1 = 2;
m2 = 2;
m3 = 2;

parameters = [g L1 L2 L3]';
mass = [m1 m2 m3]';

%% CONFIGURE SIMULATION

tf = 10;

initial_coordinates = [0.6; 0.3; -0.8] + pi;  %[th; phi; xi]
initial_velocities  = [0 0 0]';  %[dth; dphi; dxi]
init_state = [initial_coordinates; initial_velocities];

%% Simulate
Simulate_EL


%% Animate

close all;


%Pendulum 1:
obj.pend_1.type = {'line','ball'};
obj.pend_1.a = @(q) [0;0];
obj.pend_1.b = @(q) L1*arm(q(1)+pi/2);
obj.pend_1.color = GetColorCode('g',0.8);
obj.pend_1.c = obj.pend_1.b;
obj.pend_1.r = .05;

%Pendulum 2:
obj.pend_2.type = {'line','ball'};
obj.pend_2.a = @(q) obj.pend_1.b(q);
obj.pend_2.b = @(q) obj.pend_2.a(q) + L2*arm(q(2)+pi/2);
obj.pend_2.color = GetColorCode('b',0.8);
obj.pend_2.c = obj.pend_2.b;
obj.pend_2.r = .05;

%Pendulum 3:
obj.pend_3.type = {'line','ball'};
obj.pend_3.a = @(q) obj.pend_2.b(q);
obj.pend_3.b = @(q) obj.pend_3.a(q) + L3*arm(q(3)+pi/2);
obj.pend_3.color = GetColorCode('p',0.8);
obj.pend_3.c = obj.pend_3.b;
obj.pend_3.r = .05;


% CONFIGURE ANIMATION:
%Obs! the format ratio should be adjusted to your screen/figure. Usually
%one of these below will be fine.

%aspect = 5/4;
aspect = 5/4*1.55; %(or 1920/1080)
%aspect = 5/4*0.75;

% CONFIGURE ANIMATION:

% frame:
    config.framecenter = [0 -1];
    config.frameheight = 2*(L1+L2+L3); 
    config.aspect = 1920/1080;  
% figure:
    config.figureheight = 1080;
    config.figurelocation = [0 0];
% General:
    config.simspeed = 1; % slow/fast motion
    config.tf = tf;
    config.grid = 'on';
    config.enterToStart = 1;  %turn this off(to "0") to start immediately after running this sub-section
    config.video.enable = "off";

config.simspeed = 1; % slow/fast motion
config.tf = tf;
config.grid = 'on';
config.enterToStart = 1;  %turn this off(to "0") to start immediately after running this sub-section
Animate(tsim,xsim(:,1:end/2),obj,config);



