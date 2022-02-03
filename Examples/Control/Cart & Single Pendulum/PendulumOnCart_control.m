
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

%% CONFIGURE SIMULATION

tf = 20; % Duration

initial_coordinates = [0.1; 0.5];  %[th ,x]
initial_velocities  = [0;0];    %[dth,dx]
init_state = [initial_coordinates; initial_velocities];

%% CONFIGURE CONTROLELR:
%% type
controller.type = "PD"; % select PD controller
controller.FL = "off";   % use feedback linearization

%%%% Choose ONE of the configurations below:

%% Controller 1) Acutation on both cart and pendulum (PD on both terms seperately)
controller.PD.C = eye(2);   % measure both signals:  e = C*(ref - q)
controller.PD.Kp = [100 0;
                    0 10];  % u = Kp*e + Kd*de
controller.PD.Kd = [1 0;
                    0 1];
controller.B = eye(2); % Since all variables alreasy have an assiciated force in u, no need to reshape u into forces Q = B*u .
                       % one can also simply not define B, but here we wish to override it, in
                       % case you rund any of the other configurations, and come back to this one
values.umax = 10;      %Also not necessary, but is here to override any values given to this in the other controllers
                       %(when only giving ONE max value, this becomes the maximum for all control signals)

% NOTE: It is easy to stabilize fully actuated systems

%% Controller 2) Actuation only on cart (SISO PD on cart)

controller.PD.C = [0 1];    % e = C*(ref - q) (only measure x)
controller.PD.Kp = 100;      % u = Kp*e + Kd*de
controller.PD.Kd = 10;
controller.B = [0; 1];      % The control signal should be applied to x, and not th
values.umax = 10;           % Some maximum force our actuator can apply to the cart

%NOTE: This stabilizes x, but not th...

%% Controller 3) Actuation only on cart (measure th, control x)

controller.PD.C = [1 0];    % e = C*(ref - q) (only measure x)
controller.PD.Kp = 50;      % u = Kp*e + Kd*de
controller.PD.Kd = 5;
controller.B = [0; 1];      % The control signal should be applied to x, and not th
values.umax = 10;           % Some maximum force our actuator can apply to the cart

%NOTE: This stabilizes th, by only using actuation on x!

%% Controller 4) Actuation only on cart (measure both coordinate states th and x, and control x)

controller.PD.C = eye(2);
controller.PD.Kp = [50 -1];
controller.PD.Kd = [5 -1];
controller.B = [0; 1];
values.umax = 10;

% NOTE! By using the th-stabilization, and adding a weak PD on the cart
% position, we were able to stabilize this unstable cart-pendulum system(upwards),
% by only controlling the cart!!
% (see a similar 'dominant-gain' structure in the triple pendulum example)

%% SIMULATE
Simulate_EL

%% ANIMATE

close all;
clear("obj")
clear("config")

box_width = 0.4;
box_angle = pi/6;

nq = length(xsim(1,1:end/2));
nu = length(usim(1,:));

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

if length(usim(1,:)) > 1
%Input torque on Pundulum:
obj.input_pend.type = 'arrow';
obj.input_pend.tail = @(q) obj.pend.b(q);
obj.input_pend.head = @(q) obj.pend.b(q) + 0.1*q(nq+1)*arm(q(1)+pi);
obj.input_pend.color = GetColorCode('r',0.9);
end

%Input force on Cart:
obj.input_cart.type = 'arrow';
obj.input_cart.head = @(q) [q(2);0] - box_width*arm(box_angle).*[1;0]*(1 - 2*(q(nq+nu)<0));
obj.input_cart.tail = @(q) obj.input_cart.head(q) - 0.1*arm(0)*q(nq+nu);
obj.input_cart.color = GetColorCode('r',1.1);

% CONFIGURE ANIMATION:

%Adjust formatratio to your liking
formatRatio = 5/4;
%formatRatio = 5/4*1.55;
%formatRatio = 5/4*0.75;
lift = 0;
shift = 0;
height = 4*L;
width = height*formatRatio;
config.axis = [-width/2 width/2 -height/2 height/2] + [shift shift lift lift];


%Configure:
config.frameheight = 4*L; 
config.enterToStart = 1;

% Run animation
Animate(tsim,[xsim(:,1:end/2) usim],obj,config);



























