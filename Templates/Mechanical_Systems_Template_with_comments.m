
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
% controller.type = "KF";     %(The standard: u = -K*q + F*r)
% controller.type = "EXT";    %(if you want to use an externally defined controller)

%%% To use Feedback LInearization:
%controller.FL = "on"

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
controller.KF.K = place(double(S.LA),double(S.LB),poles);

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

%% Animate (some example objects)
clear("obj"); clear("config"); close all;

%Coordinates used to define positions  in the animation:
AnimationValues = [xsim(:,1:end/2) usim' Ref];


%%% DEFINE OBJECTS TO ANIMATE:
%format:
%{
collect all objects to animate into struct "obj"

An object looks like:
obj.ObjectName.type = ' ';
obj.ObjectName.*property_1* = @() ...
obj.ObjectName.*property_2* = @() ...
obj.ObjectName.*property_3* =     ...
(not required):
obj.ObjectName.color = [*RGB*] or GetColorCode(' blue/red/green/... ',brightness)
obj.ObjectName.disable = 0 / 1;

All object are also available as "timed_*type*", where all function handles
take time as argument.
%}


%Box (cart):
obj.myBox.type = 'box';
obj.myBox.def = "CD";  %Define box with a center and length of diagonal(can also use: Center/Any corner(CC), Two Opposite Corners(2C) )
obj.myBox.B1 = @(q) [ ; ]; %center
obj.myBox.B2 = @(q) ;      %diagonal
obj.myBox.color = GetColorCode('p',1.1); %give color to box
obj.myBox.disable = ; % 1 - diables box(is not animated)

%Line:
obj.myLine.type = 'line';
obj.myLine.a = @(q) [ ; ]; % start of line
obj.myLine.b = @(q) [ ; ]; % end of line
obj.myLine.color = GetColorCode('i'); %select a color
obj.myLine.disable = ;

%Ball (circle):
obj.myBall.type = 'ball';
obj.myBall.c = [ ; ]; %center
obj.myBall.r = .05;   %radius
obj.myBall.color = GetColorCode('blue',0.8);
obj.myBall.disable = ;

%Arrow (vectors):
obj.myArrow.type = 'arrow';
obj.myArrow.tail = @(q) [ ; ];
obj.myArrow.head = @(q) [ ; ];
obj.myArrow.color = GetColorCode('r',0.9);

%Point (dot):
obj.myPoint.type = 'point';
obj.myPoint.p = @(q) [ ; ]; %placement
obj.myPoint.color = [0 0 0]; % black is default (can also use GetColorCode('k/black'))

%Circular Arrow (moment/torque):
obj.myCircleArrow.type = 'moment';
obj.myCircleArrow.target = @(q) [ ; ];    %center of circle
obj.myCircleArrow.magnitude = @(q) [ ; ]; %radius
obj.myCircleArrow.minr = ;                %smallest radius(radius will be equal to: minr + magnitude)
obj.myCircleArrow.color = GetColorCode('r',0.8); %Select main color
obj.myCircleArrow.color2 = GetColorCode('b',1.2);%Color when going in reversed direction
obj.myCircleArrow.thickness = @(q)  ; %LineWidth

%%%% Combined objects:

%Pendulum:
obj.myPendulum.type = {'line','ball'};
obj.myPendulum.a = @(q) [ ; ];
obj.myPendulum.b = @(q) obj.ref_1.a(q) + [ ; ];
obj.myPendulum.c = obj.ref_1.b;
obj.myPendulum.r = ;
obj.myPendulum.color = GetColorCode('o');
obj.myPendulum.disable = ;



%%%%%%%%%% CONFIGURE ANIMATION:

%the only necessary configuration to give is either axis or frameheight, the rest is optional.

%%%%%% frame:
% alternative 1)  (prioritized)
    config.axis = [ , , , ];    %(when defined, this takes priority)
% alternative 2)
    config.framecenter = [ , ];
    config.frameheight = ; 
    % and either of: (
    config.aspect = 1920/1080;  
    config.framewitdth = ; %(prioritized when both aspect and )
    %)

%%%%%% figure:
    %either of:
    config.position = [1 1 1920 1080]; %(this takes priority)
    %
    config.figureheight = 1080;  % gets same aspect as frame
    config.figurelocation = [0 0];
    %

%other configuration options
config.simspeed = 1;  % set the speed of the animation relative to real time(2 - twice as fast motion, 0.5 - half the speed(slow motion))
config.tf = tf;       % set the time to stop animation(default: animates the enitre timeseries "tsim")
config.grid = 'on';   % same "grid" property as normal plots
config.enterToStart = 1; % set this to make the animation wait for you to press enter before stating

% Animate:
Animate(tsim,AnimationValues,obj,config);



%% Save Animation as video:

config.video.enable = "on";     %make "Animation" produce a video instead of animating
config.Video.resolution = 0.01; %Set the number of seconds that passes between frames(seconds in the animation, that is)

config.video.profile = 'MPEG-4';        % some video prfile (from "VideoWriter" by matlab)
config.video.LosslessCompression = 1;   % toggle on or off(only applies to some dataformats)
config.video.CompressionRatio = 2;      % set compression ratio(only applies to some dataformats)
  
Animate(tsim,[xsim(:,1:end/2) usim],obj,config);  %Use "Animate" to save a video

% The video appears in the current directory, and is called:
% "AnimationVideo". Change this name if you want to save it, so that is is
% not overridden by the next video you make.



