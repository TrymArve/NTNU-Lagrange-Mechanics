
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


%% 1) EXAMPLE 1: Only control th and phi (only have motors on the two first joints)
clear("controller")

tf = 20;

initial_coordinates = [0.1; 0.1; -0.1];  %[th; phi; xi]
initial_velocities  = [0 0 0]';  %[dth; dphi; dxi]
init_state = [initial_coordinates; initial_velocities];

controller.type = "PD";

controller.PD.C = eye(3);
controller.PD.Kp = [50 -35   0;
                     0  100 -1000];
controller.PD.Kd = [-10 -10   0;
                     0  -10  -100];

controller.B = [1 0; 0 1; 0 0];
values.umax = 1000;


%% 2) EXAMPLE 2: Only control th! (only have motor on first joint)
clear("controller")

tf = 20;
%{ 
Other stable initial coordinates:
initial_coordinates = [0.1; 0.1; -0.01];    %[th; phi; xi]
initial_coordinates = [0.2; -0.1; -0.01];    %[th; phi; xi]
initial_coordinates = [0.4; 0; -0.02];    %[th; phi; xi]
initial_coordinates = [0.4; -0.05; -0.02];    %[th; phi; xi]
%}
initial_coordinates = [0.1; -0.1; -0.02];    %[th; phi; xi]

initial_velocities  = [0 0 0]';             %[dth; dphi; dxi]
init_state = [initial_coordinates; initial_velocities];

controller.type = "PD";

controller.PD.C = eye(3);
controller.PD.Kp = [130  -1000  2000];
controller.PD.Kd = [350  -100  1800];

controller.B = [1 0 0]';
values.umax = 10000;


%% No Controller:
clear("controller")
controller.type = "off";

%% Simulate
Simulate_EL


%% Animate
clc;close all;clear("obj");clear("config");

nq = length(xsim(1,1:end/2));

radius = 0.05;
magnitudeScaleing = 0.2;

%Pendulum 1:
obj.pend_1.type = {'line','ball'};
obj.pend_1.a = @(q) [0;0];
obj.pend_1.b = @(q) L1*arm(q(1)+pi/2);
obj.pend_1.color = GetColorCode('g',0.8);
obj.pend_1.c = obj.pend_1.b;
obj.pend_1.r = radius;

%Pendulum 2:
obj.pend_2.type = {'line','ball'};
obj.pend_2.a = @(q) obj.pend_1.b(q);
obj.pend_2.b = @(q) obj.pend_2.a(q) + L2*arm(q(2)+pi/2);
obj.pend_2.color = GetColorCode('b',0.8);
obj.pend_2.c = obj.pend_2.b;
obj.pend_2.r = radius;

%Pendulum 3:
obj.pend_3.type = {'line','ball'};
obj.pend_3.a = @(q) obj.pend_2.b(q);
obj.pend_3.b = @(q) obj.pend_3.a(q) + L3*arm(q(3)+pi/2);
obj.pend_3.color = GetColorCode('p');
obj.pend_3.c = obj.pend_3.b;
obj.pend_3.r = radius;

% Moment 1:
obj.input_1.type = 'moment';
obj.input_1.target = @(q) obj.pend_1.a(q);
obj.input_1.magnitude = @(q) 0.1*magnitudeScaleing*q(nq+1);
obj.input_1.minr = 0.01;
obj.input_1.color = GetColorCode('r',0.8);
obj.input_1.color2 = GetColorCode('b',1.2);
obj.input_1.thickness = @(q)  0.01 + abs(magnitudeScaleing*0.1*q(nq+1));

% Moment 2:
if length(usim(1,:)) > 1
obj.input_2.type = 'moment';
obj.input_2.target = @(q) obj.pend_2.a(q);
obj.input_2.magnitude = @(q) 0.1*magnitudeScaleing*q(nq+2);
obj.input_2.minr = 0.01;
obj.input_2.color = GetColorCode('r',0.8);
obj.input_2.color2 = GetColorCode('b',1.2);
obj.input_2.thickness = @(q)  0.01 + abs(magnitudeScaleing*0.1*q(nq+1));
end



% CONFIGURE ANIMATION:

%%%%%% frame:
% alternative 1)  (prioritized)
    config.framecenter = [0 1];
    %config.axis = [1 10 2 20];  
% alternative 2)
    config.framecenter = [0 1];
    config.frameheight = 2*(L1+L2+L3); 
    % and either of: (
    config.aspect = 1920/1080; %(or 5/4 or 5/4*0.75 or 5/4*1.55, change this to make the frame not skewed) 
    config.framewitdth = 2*(L1+L2+L3); %(priority)
    %)

%%%%%% figure:
    %either of:
    config.position = [1 1 1920 1080];
    %
    config.figureheight = 1080;
    config.figurelocation = [0 0];
    %

%%%%% General:
config.simspeed = 1; % slow/fast motion
config.tf = tf;
config.grid = 'on';
config.enterToStart = 1;  %turn this off(to "0") to start immediately after running this sub-section
config.video.enable = "off";

%%%%% ANIMATE:
Animate(tsim,[xsim(:,1:end/2) usim],obj,config);

%% Save Animation as video

config.video.enable = "on";
config.video.profile = 'MPEG-4';
config.video.LosslessCompression = 1;
config.video.CompressionRatio = 2;
config.Video.resolution = 0.01;
Animate(tsim,[xsim(:,1:end/2) usim],obj,config);

















