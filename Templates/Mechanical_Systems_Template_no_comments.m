
%% Define System (2D)
clc;clear;close all;


S = MakeSymbolicGenCoords({'var1','var2','var3','...'},{'par1','par2','par3','...'}); 

pos_1   = [ , ];
pos_2 = [ , ];
S.pos = {pos_1,pos_2};

S = MakeLagrange(S);

%% Give values to the various system parameters and masses:

% Parameters:
g = 9.81;
parameter_1 = ;
parameter_2 = ;
parameters = [g;parameter_1;parameter_2 ];

% Masses:
m1 = ;
m2 = ;
masses = [m1;m2];


%% Configure the Simulation:

tf = ;

initial_coordinates = [ ; ; ; ];
initial_velocities  = [ ; ; ; ];
init_state = [initial_coordinates; initial_velocities];

%% Configure Controller:

% controller.type = "off";
% controller.type = "PD";  
% controller.type = "FLPD"; 
% controller.type = "KF";    
% controller.type = "FLKF"; 
% controller.type = "EXT";   

controller.B =  ;

%% Design the controller you want to use:
%% PD:
controller.PD.C  = [ , , , ]; 
controller.PD.Kp = ;
controller.PD.Kd = ; 

%% KF:  (always linearizes about reference point)

%%%%%%%%% Get Symbolic Linearization:
S = LinearizeEL(S);
LA = @(q,dq) LinA(q,dq,mass,parameters); 
LB = @(q,dq) LinB(q,dq,mass,parameters);

poles = [ , , , ];
controller.KF.K = @(t,q,dq) place(LA(q,dq),LB(q,dq),poles);

%% KF: Linearize about fixed point:

%%%%%%%%% Get Numeric Linearization (about the point: [q0,dq0]):
q0  = [ ; ; ; ];
dq0 = [ ; ; ; ];
poles = [ , , , ];
S = LinearizeEL(S,q0,dq0,masses,parameters);
controller.KF.K = place(S.LA,S.LB,poles);

%% EXT:

controller.EXT.source = "timeseries";
controller.EXT.timeseries = [ ... ; ...];

controller.EXT.source = "function";
controller.EXT.function = @(t,q,dq) ... ;



%% Set the Reference Trajectory:
freq = 1;
amp = 1;
controller.ref   = @(t) [ amp*sin(freq.*t);        0.*t; 0.*t; 0.*t];
controller.dref  = @(t) [ amp*freq*cos(freq.*t);   0.*t; 0.*t; 0.*t];
controller.ddref = @(t) [-amp*freq^2*sin(freq.*t); 0.*t; 0.*t; 0.*t];

%% Some values the integrator needs

% VALUES:
values.umax = [ ; ; ; ];


%% SIMULATION  
Simulate_EL;

%% Animate  (Various Object Examples Included)

close all;
clear("obj")

AnimationValues = [xsim(:,1:end/2) usim' Ref];
startU   = size(xsim(1,1:end/2),2);
startRef = startU + size(usim',2);

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

%Input torque on Pundulum 1:
obj.input_pend_1.type = 'arrow';
obj.input_pend_1.tail = @(q) obj.pend_1.b(q);
obj.input_pend_1.head = @(q) obj.pend_1.b(q) + 0.1*q(startU+2)*arm(q(2)+pi);
obj.input_pend_1.color = GetColorCode('r',0.9);

%Input torque on Pundulum 2:
obj.input_pend_2.type = 'arrow';
obj.input_pend_2.tail = @(q) obj.pend_2.b(q);
obj.input_pend_2.head = @(q) obj.pend_2.b(q) + 0.1*q(startU+3)*arm(q(3)+pi);
obj.input_pend_2.color = GetColorCode('r',0.9);

%Input torque on Pundulum 2:
obj.input_pend_3.type = 'arrow';
obj.input_pend_3.tail = @(q) obj.pend_3.b(q);
obj.input_pend_3.head = @(q) obj.pend_3.b(q) + 0.1*q(startU+4)*arm(q(4)+pi);
obj.input_pend_3.color = GetColorCode('r',0.9);


% CONFIGURE ANIMATION:
formatRatio = 5/4;
%formatRatio = 5/4*1.55;
%formatRatio = 5/4*0.75;
lift = -.5;
shift = 0;
height = 3.5*(L1+L2);
width = height*formatRatio;
config.axis = [-width/2 width/2 -height/2 height/2] + [shift shift lift lift];
config.simspeed = 1;
config.tf = tf;
config.grid = 'on';
config.enterToStart = 1;
Animate(tsim,AnimationValues,obj,config);



