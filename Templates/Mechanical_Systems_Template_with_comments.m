
%% Define System (2D)
clc;clear;close all;

% 1) Define Generalized Coordinates Variables and System Parameters:
S = MakeSymbolicGenCoords({'var1','var2','var3','...'},{'par1','par2','par3','...'}); 
% Creates a struct containing fields: 'q', 'dq', 'ddq' which are vectors
% containing the symbolic variables and their single and double derivatives
% (also a vector containing the coordinates as symbolic functions of time: q(t))
% This struct will be expanded to contain your system


% 2) Define Positions of point masses:
pos_1   = [ , ];
pos_2 = [ , ];
pos_3 = [ , ];
pos_4 = [ , ];
pos_ ... ;
S.pos = {pos_1,pos_2,pos_3pos_4, ... }; % add a cell-array containing the positions to the struct
% Respective masses are automatically defined as [m1;m2;m3;...] as symbolic vaiables



% 3) Compute the Euler-Lagrange equation:
S = MakeLagrange(S);
% (The idea is to contain all system properties in this struct)
% (use: "S.help" to get a description of fields)


%% Give values to the various system parameters and masses:

% Parameters:
g = 9.81; %(Required) Gravity is automatically the first element of the parameters vector.
parameter_1 = ;
parameter_2 = ;
parameter_3 = ;
parameter_ ... ;
parameters = [g;parameter_1;parameter_2;parameter_3, ... ]; %Put g first, then the rest in the same order as you put them into "MakeSymbolicGenCoords"

% Masses:
m1 = ;                 %The mass at the fisrt  position desribed above: pos_1
m2 = ;               %The mass at the second position desribed above: pos_2
m3 = ;                 %The mass at the third  position desribed above: pos_3
m4 = ;               %The mass at the fourth position desribed above: pos_4
m ... ;
masses = [m1;m2;m3;m4; ... ];


%% Configure the Simulation:

%%% Time of Finish: (for how long to simulate in seconds)
tf = ;  %Put the duration of the simulation in seconds

%%% Initial State for the simulation:
initial_coordinates = [ ; ; ; ];
initial_velocities  = [ ; ; ; ];
init_state = [initial_coordinates; initial_velocities];

%% Configure Controller:

%%% Select what controller to use:  (uncomment the one you want to use)
% controller.type = "off";    %(no controller used)
% controller.type = "PD";     %(PD-controller. Meant for SISO systems)
% controller.type = "FLPD";   %(Feedback Linearization with PD. Meant for SISO)
% controller.type = "KF";     %(The standard: u = -K*q + F*r)
% controller.type = "FLKF";   %(Feedback Linearization with "KF")
% controller.type = "EXT";    %(if you want to use an externally defined controller)

controller.B =  ; % for transforming the inputsignal into the force on the system. 
% i.e. M*ddq+C*dq+G = Q = B*u
% Q is the force on the system
% u is the input vector 
% (if you don't have actuation available in the second link of an arm, but do have actuation in the first and third, the integrator still needs to know that the force is zero in the second link:
% u = [link1; link3] ---> Q = [link1; 0; link3] 
% then:  B = [1   0;
%             0   0;
%             0   1] )

%% Design the controller you want to use:
%% PD:
controller.PD.C  = [ , , , ]; % Measure output (y = C*q) --> (e = C*(q - ref))
controller.PD.Kp = ; %Proportional gain
controller.PD.Kd = ; %Derivative gain

%% KF:  (always linearizes about reference point)

%%%%%%%%% Get Symbolic Linearization:
S = LinearizeEL(S); %Provides symbolic LA and LB used in the linearization:   [dq;ddq] = LA*[q;dq] + LB*u
LA = @(q,dq) LinA(q,dq,mass,parameters); % Obtain Linearization for point [q;dq]
LB = @(q,dq) LinB(q,dq,mass,parameters);

poles = [ , , , ];
controller.KF.K = @(t,q,dq) place(LA(q,dq),LB(q,dq),poles);
% We don't need a feedforward gain F, since we're using dynamic linearization

%% KF: Linearize about fixed point:

%%%%%%%%% Get Numeric Linearization (about the point: [q0,dq0]):
q0  = [ ; ; ; ];
dq0 = [ ; ; ; ];
poles = [ , , , ];
S = LinearizeEL(S,q0,dq0,masses,parameters);
controller.KF.K = place(S.LA,S.LB,poles);

%% EXT:
% Select controller source:

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
 % Add ".*t" so that all elements are of the same size when inputing a timeseries

%% Some values the integrator needs

% VALUES:
values.umax = [ ; ; ; ];    % Maximum allowed vaules on control signals (before transformation by B) (not reqiured)



%% SIMULATION  

%Simply run this to simulate the system:
Simulate_EL;

%% Animate

close all;
scale = 1;


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
formatRatio = 5/4*1.55;
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



