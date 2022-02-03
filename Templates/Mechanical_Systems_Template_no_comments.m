
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

%% Animate (some example objects)
clear("obj"); clear("config"); close all;

AnimationValues = [xsim(:,1:end/2) usim' Ref];

%Box (cart):
obj.myBox.type = 'box';
obj.myBox.def = "CD";
obj.myBox.B1 = @(q) [ ; ];
obj.myBox.B2 = @(q) ;     
obj.myBox.color = GetColorCode('p',1.1);
obj.myBox.disable = ; 

%Line:
obj.myLine.type = 'line';
obj.myLine.a = @(q) [ ; ];
obj.myLine.b = @(q) [ ; ]; 
obj.myLine.color = GetColorCode('i');
obj.myLine.disable = ;

%Ball (circle):
obj.myBall.type = 'ball';
obj.myBall.c = [ ; ];
obj.myBall.r = .05;   
obj.myBall.color = GetColorCode('blue',0.8);
obj.myBall.disable = ;

%Arrow (vectors):
obj.myArrow.type = 'arrow';
obj.myArrow.tail = @(q) [ ; ];
obj.myArrow.head = @(q) [ ; ];
obj.myArrow.color = GetColorCode('r',0.9);

%Point (dot):
obj.myPoint.type = 'point';
obj.myPoint.p = @(q) [ ; ]; 
obj.myPoint.color = [0 0 0]; 

%Circular Arrow (moment/torque):
obj.myCircleArrow.type = 'moment';
obj.myCircleArrow.target = @(q) [ ; ];  
obj.myCircleArrow.magnitude = @(q) [ ; ];
obj.myCircleArrow.minr = ;             
obj.myCircleArrow.color = GetColorCode('r',0.8); 
obj.myCircleArrow.color2 = GetColorCode('b',1.2);
obj.myCircleArrow.thickness = @(q)  ; 

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

    config.axis = [ , , , ];
    config.framecenter = [ , ];
    config.frameheight = ; 
    config.aspect = 1920/1080;  
    config.framewitdth = ;

    config.position = [1 1 1920 1080]; %(this takes priority)
    config.figureheight = 1080;  % gets same aspect as frame
    config.figurelocation = [0 0];

config.simspeed = 1; 
config.tf = tf;     
config.grid = 'on';  
config.enterToStart = 1; 

% Animate:
Animate(tsim,AnimationValues,obj,config);

%% Save Animation as video:

config.video.enable = "on";    
config.Video.resolution = 0.01;

config.video.profile = 'MPEG-4';      
config.video.LosslessCompression = 1; 
config.video.CompressionRatio = 2;
  
Animate(tsim,[xsim(:,1:end/2) usim],obj,config);

