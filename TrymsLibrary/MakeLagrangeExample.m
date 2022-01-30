% MAKE-LAGRANGE EXAMPLE

clc;clear;close all;




%%%%%%% REQUIRED:

% Generalized Coordinates
GenCoords = {'x1', 'x2', 'e1', 'e2'};
mysystem  = MakeSymbolicGenCoords(GenCoords);

syms g  % is automatically added if not included here
syms L  % other relevant parameters used in eucledean positions
mysystem.par = [g L];

% Define Euclidean positions based on the generalized coordinates(functions of time):
%(ex: inverted pendulum on a puck)
P_mount = [x1; x2; 0];
P_pend = P_mount + L*[cos(e1)*sin(e2);
                      sin(e1)*sin(e2);
                      cos(e2)];
mysystem.pos = {P_mount, P_pend}; % add all positions here        
           




%%%%%%% OPTIONAL: (If not defined, the function behaves as if they are 0)

% mysystem.subDiffs = 1; %Provide the sub-differentiations with the output-struct(dL_q, dL_dq, dL_dq_dt)
                       %This also prodives energies og specific masses
 mysystem.extractMCG  = 1; %Extract the matrices M,C,G.  ("sytem model": M*ddq + C*dq + G = 0)
 mysystem.compareMCG  = 1; %Compares them to the Lagrange equations(to make sure they are correct)
 mysystem.nofuncMCG   = 1; %Does NOT create a matlab function to find numeric values for M,C,G for a given (q,dq)
% mysystem.extractW = 1; %Extracts momentum-matrix (note: W is the same as M)
% mysystem.nofuncW  = 1; %Makes a matlabFunction to return the
%                         momentum-matrix and the RHS equations
%%%% Masses of the system: (is these are not defined, they are generated automatically)
% syms M m   %Should have same order ar the positions(-->pos(1) describes pos of mass(1))
% mysystem.mass = {M ,m}; 
% (NOTE: if no masses are defined here, the variables "m1","m2",...,"mN" are reserved)





%%%%%%% MAKE THE LAGRANGE SYSTEM:
Lag = MakeLagrange(mysystem);

MakeLagrange('help')    % Displays desciptions about configuration variables
Lag.help                % Displays some descriptive text about the variables
Lag                     % Displays the available variables

%%%%%%% FORCES ARE USED AFTER THIS FUNCTION AND SHOULD NOT BE INPUTS
% Forces (also consistent indecies)
% syms u1 u2
% Q = {u1, u2, 0,  0};  %Generalized forces

%% OTHER EXAMPLES:

%% System from Homework 11
clc;clear;close all;
% Based on "MakeLagrangeExample"


% Generalized Coordinates
GenCoords = {'x', 'th'};
mysystem  = MakeSymbolicGenCoords(GenCoords);

syms L  % other relevant parameters used in eucledean positions
mysystem.par = L;

% Define Euclidean positions based on the generalized coordinates(functions of time):
P_box = [x; 0];
P_pend = P_box + L*[sin(th);
                    -cos(th)];
mysystem.pos = {P_box, P_pend}; % add all positions here        
           
%config
mysystem.extractW = 1;
mysystem.extractMCG = 1;

Lag = MakeLagrange(mysystem)
Lag.EL
Lag.W
Lag.M
Lag.C
Lag.G

%% System from Homework 9 (Modeling Inertia!!)
clc;clear;close all;

GenCoords = {'q1','q2'};
S = MakeSymbolicGenCoords(GenCoords);
syms Lc1 L1
S.par = [Lc1,L1];

pos_rod = Lc1*[cos(q1);sin(q1)];
pos_rod_1 = pos_rod + [cos(q1);sin(q1)];
pos_rod_2 = pos_rod - [cos(q1);sin(q1)];
pos_wheel = L1*[cos(q1);sin(q1)];
pos_wheel_1 = pos_wheel + 1*[cos(q1+q2); sin(q1+q2)];
pos_wheel_2 = pos_wheel - 1*[cos(q1+q2); sin(q1+q2)];
S.pos = {pos_rod,pos_rod_1,pos_rod_2,pos_wheel,pos_wheel_1,pos_wheel_2};

S.extractMCG = 1;
S.extractW = 1;
S = MakeLagrange(S);
S.help;

EL = RepairInertia(S.EL, {'m1','m2','m3'},{'I1'});
EL = RepairInertia(EL, {'m4','m5','m6'},{'I2','m2'});

%% Help

MakeLagrange('help')
MakeLagrange('config')
MakeLagrange('fields')
MakeLagrange('info')
MakeLagrange('funcs')

%% EXTRACTOR TEST
clc;clear;
syms x y z dx dy dz ddx ddy ddz
EL = [dx^2 + dy^2 + dz^2 + dx*dy + dy*z + dx*dz ;
      4*dx^2 + 5 + dy + 8*dz;
      sin(dx) + dz ]
[SM,M,C,G] = ExtractMCG(EL,{'x','y','z'});
