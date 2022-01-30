function[Lagrange] = MakeLagrange(sys)
if ~isstruct(sys)
    if string(sys) == "help"
        disp('Try: "config" - for configuration options');
        disp('Or:  "fields" - for the required fields for making the Lagrange Mechanics');
        disp('Or:  "info"  - for general infomation')
        disp('Or:  "funcs" - for related functions')
        return;
    elseif string(sys) == "config"
        config.subDiffs   = "The various partial differentiations of the Lagrangian are provided in the output-struct (dL_q, dL_dq, dL_dq_dt) You also get energies of individual masses";
        config.extractMCG = "The system model matrices on the form (M(q)*ddq + C(q,dq)*dq + G(q) = u) will be extracted from the Euler Lagrange equations and given as output";
        config.nofuncMCG  = "The model matrices are NOT made available via a function MCG(q,dq,mass,par). (such function is made by default when extractMCG is enabled)";
        config.extractW   = "The the momentum-matrix on the form (M(q)*ddq + C(q,dq)*dq + G(q) = u) will be extracted from the Euler Lagrange equations and given as output";
        config.nofuncW    = "The momentum-matrix and RHS is NOT made available via a function W_RHS(q,dq,mass,par). (such function is made by default)";
        
        config
        return;
    elseif string(sys) == 'fields'
        fields.Required = "REQUIRED:";
        fields.GenCoords = "Cell array of character-strings of the variable names"
        fields.q = "Symbolic variables of the generalized coordinates";
        fields.dq = "Symbolic variables of the time derivatives of q";
        fields.ddq = "Symbolic variables of the double time derivatives of q";
        fields.dt = "Generalized coordinates q, but as symbolic functions of time. ex: 'syms x1(t)'";
        fields.pos = "A cell array of the Eucledean position vectors of each mass in the system";
        fields.par = "An array of symbolic parameters used the position vectors. ex: lengths of arms";
        fields.hint = "Use the function 'MakeSymbolicGenCoords({'q1','q2',..})' to create the symbolic variables needed";

        fields.O = "------------";
        fields.Optional = "OPTIONAL:";
        fields.mass = "Assign your own names to the mass-variables. use ex: Sys.mass = {M,m,m2,M5};(cells containing SYMBOLIC VARIABLES). indices correspond to the positions";
        fields.q = "(alternative to GenCoords)Provide the necessary symbolic coordinate variables istead of GenCoords(see info in function)";
        fields.MakeGenCoords = "(if alternative coordinates)Use 'MakeSymbolicGenCoords()' to generate the required coordinates";
        
        fields
        return;
    elseif string(sys) == 'info'
        info.W_M = "The matrix W and M are the same matrix, but are calucalted differently and are inspired by different approaches";
        info.W_RHS = "The W and RHS system is intended for easily simlulating the system";
        info.MCG = "The MCG system is intended to help design model based controllers";
        info.support = "There are a couple of functions created to help you use this function more easily(see 'funcs')";
        info.inertia = "This function can ONLY model point masses, HOWEVER!, you can easily add inertia by adding two other masses in opposite directions from the center of mass of the body. (such that there are three masses per body) Then 'Repair' the model afterwards with 'RepairInertia()' ";
        info.g       = "OBS! gravity must always point downwards in the LAST euclidean dimension. (if to dimentional: [0, -g] , three: [0,0,-g], four: [0,0,0,-g], etc.)";

        
        info
        return;
    elseif string(sys) == 'funcs'
        funcs.MakeSymbolicGenCoords = "Makes the necessary symbolic variables based on a list of names. outputs a struct. Add this struct to the system struct(merge, not substruct)";
        funcs.mergestructs = "Use this to merge structs";
        funcs.RepaiInertia = "Use this to swap out 'side-masses' with an inertia variable, to obtain the proper model";
        funcs.ExtractMCG = "Extracts the symbolic M,C,G matrices of the input equations of motions.(for use after RepairInertia)";
        
        funcs
        return;
    end
end


% Use "Lagrange.help" for info on the various variables

% REQUIRED INPUTS:
% sys.par  - vector of symbolic parameters needed in the convertion from
%            generalized coordinates to Eucledean coordinates
%            ("g" is reseved for gravity, and is added automatically if it is not already included)
% sys.GenCoords - Cell array of character-strings of the variable names
% sys.qt   - cell array of symbolic generalized coordinates (MUST be functions of time: "x1(t)") 
%            (the ones used in defining euclidean positions)
% sys.q    - same as qt, but NOT functions of time...
% sys.dq   - cell array of symbolic derivatives (not functions of time)
% sys.ddq  - cell array of symbolic double derivatives (not functions of time)
% sys.pos  - cell array of euclidean positions of the listed point masses (built with q-values)

% OPTIONAL INPUTS: (if not defined, they default as 0)
% sys.mass - cell array of masses of corresponding position
% sys.subDiffs    - if == 1: The various partial differentiations of the
%                            Lagrangian are provided in the output-struct
%                            (dL_q, dL_dq, dL_dq_dt)
%                            You also get energies of individual masses

% sys.extractMCG  - if == 1: The system model matrices on the form (M(q)*ddq + C(q,dq)*dq + G(q) = u)
%                            will be extracted from the Euler Lagrange
%                            equations and given as output

% sys.compareMCG  - if == 1: The extracted model will be compared to the
%                            EL-equations and the result is printed as a
%                            message. (to verify the validity of the matrices)
%                           (if they do not match, then the system could contain 
%                            higher exponents than 1 and 2 for the ddq and dq)

% sys.nofuncMCG   - if == 1: The model matrices are NOT made available via a
%                            function "MCG(q,dq,mass,par)". 
%                            (such function is made by default)

% sys.extractW    - if == 1: The the momentum-matrix on the form (M(q)*ddq + C(q,dq)*dq + G(q) = u)
%                            will be extracted from the Euler Lagrange
%                            equations and given as output

% sys.nofuncW     - if == 1: The momentum-matrix is NOT made available via a
%                            function "W(q,dq,mass,par)". 
%                            (such function is made by default)








%%%%%%%%% SET UP SYSTEM:
syms t                 % make 'time' variable
Lagrange = sys;        % pass the input-struct back out
pos = sys.pos;
GenCoords = sys.GenCoords;
N = length(pos);  %Number of masses included in system

% Miscellaneous parameters: (what ever is used in positions)
% (ADD g to parameters if it is not there already)
g_ind = find(Lagrange.par == 'g');
if isempty(g_ind)
    syms g                           % must use g to calculate energies
    Lagrange.par = [g; Lagrange.par]; % add g to parameters
    par = Lagrange.par;
else
    par = Lagrange.par;
    g = par(g_ind);   
end
Lagrange.help.par = "Parameters needed to convert the GenCoords to Euclidean coordinates. Such as lengths of arms.";

% Masses (used in energy calculations)
if isfield(sys, 'mass')
    mass = sys.mass;
else
    m = sym('m',[1,N]);        % Make symbolic variable for each mass in the system
    for ii = 1:N
    Lagrange.mass{ii} = m(ii);
    mass = Lagrange.mass;
    end
end
Lagrange.help.mass = "Masses in the system. The corresponding masses and posistions have the same index.";

% States (required)
qt  = sys.qt;
q   = sys.q;
dq  = sys.dq;
ddq = sys.ddq;
%(Include in output)
Lag.q   = q;
Lag.dq  = dq;
Lag.ddq = ddq;


% ADD CONSTRAINTS:
if isfield(sys,'constraints')
    constraints = sys.constraints;
    lambda = sys.lambda;
else
    constraints = 0;
    lambda = 0;
end


%%%%%%%%% AUTOMATIC GENERATION OF MECHANICS:
n = length(qt);        %DoF in GenCoords
neuc = length(pos{1}); %DoF in Eucledean coordinates

% DEFINE EUCLEDEAN POSISTIONS:
pos = sys.pos;
N = length(pos);  %Number of masses included in system

for k = 1:N
    for i = 1:n
        pos{k} = subs(pos{k}, q{i}, qt{i});
        
    end
end

for k = 1:length(constraints)
    for i = 1:n
        constraints(k) = subs(constraints(k), q{i}, qt{i});
        
    end
end


% ENERGIES:
up = zeros(1,neuc); up(end) = 1; % Define the last direction as up(opposing gravity)
V = cell(N,1); P = cell(N,1); K = cell(N,1);
Ktot = 0; Ptot = 0;
for i = 1:N
V{i} = diff(pos{i}, t);             % Find velocities
P{i} = g*mass{i}*up*pos{i};         % Find potential energies
K{i} = (1/2)*mass{i}*(V{i}.'*V{i}); % Find kinetic energies

Ptot = Ptot + P{i};
Ktot = Ktot + K{i};
end

% LAGRANGIAN:
Lag = Ktot - Ptot + lambda.'*constraints;


    

% EULER-LAGRANGE EQUATIONS:
dL_q     = cell(n,1);
dL_dq    = cell(n,1);
dL_dq_dt = cell(n,1);
EulLag   = cell(n,1); EL = [];
for i = 1:n
% Differentiate wrt q:
dL_q{i} = simplify( diff(Lag, qt{i}) );

    
% Differentiate wrt dq:
dL_dq{i}    = simplify( diff(Lag     , diff(qt{i}, t)) );
dL_dq_dt{i} = simplify( diff(dL_dq{i}, t            ) );

% Euler Lagrange equation:
EulLag{i} = dL_dq_dt{i} - dL_q{i};

% Introduce derivative variables:
EL = [EL; EulLag{i}];
end


% Substitute coordinates and derivatives:

for i = 1:n
          Ptot  = simplify( subs(Ptot, diff(qt{i}, t), dq{i}) );
          Ptot  = simplify( subs(Ptot, qt{i}, q{i}) );
          Ktot  = simplify( subs(Ktot, diff(qt{i}, t), dq{i}) );
          Ktot  = simplify( subs(Ktot, qt{i}, q{i}) );

    for ii = N
          P{ii}  = simplify( subs(P{ii}.', diff(qt{i}, t), dq{i}) );
          P{ii}  = simplify( subs(P{ii}.', qt{i}, q{i}) );
          K{ii}  = simplify( subs(K{ii}.', diff(qt{i}, t), dq{i}) );
          K{ii}  = simplify( subs(K{ii}.', qt{i}, q{i}) ); 
    end
    
    if isfield(sys, 'subDiffs') && sys.subDiffs == 1
        for ii = 1:n
            dL_q{ii} = simplify(subs(dL_q{ii}, diff(qt{i}, t, t), ddq{i}));
            dL_q{ii} = simplify(subs(dL_q{ii}, diff(qt{i}, t), dq{i}));
            dL_q{ii} = simplify(subs(dL_q{ii}, qt{i}, dq{i}));

            dL_dq{ii} = simplify(subs(dL_dq{ii}, diff(qt{i}, t, t), ddq{i}));
            dL_dq{ii} = simplify(subs(dL_dq{ii}, diff(qt{i}, t), dq{i}));
            dL_dq{ii} = simplify(subs(dL_dq{ii}, qt{i}, dq{i}));

            dL_dq_dt{ii} = simplify(subs(dL_dq_dt{ii}, diff(qt{i}, t, t), ddq{i}));
            dL_dq_dt{ii} = simplify(subs(dL_dq_dt{ii}, diff(qt{i}, t), dq{i}));
            dL_dq_dt{ii} = simplify(subs(dL_dq_dt{ii}, qt{i}, dq{i}));
        end
    end
    EL = simplify(subs(EL, diff(qt{i}, t, t), ddq{i}));
    EL = simplify(subs(EL, diff(qt{i}, t), dq{i}));
    EL = simplify(subs(EL, qt{i}, q{i}));
end
    Lagrange.Pi = P;
    Lagrange.Ki = K;
    Lagrange.P  = Ptot;
    Lagrange.K  = Ktot;
    Lag = simplify(Ktot - Ptot + lambda.'*constraints);
    Lagrange.Lag = Lag;
Lagrange.help.Pi  = "The potential energies of the individual masses. Corresponding indecies to the masses in -mass-";
Lagrange.help.Ki  = "The kinetic energies of the individual masses. Corresponding indecies to the masses in -mass-";
Lagrange.help.P   = "The total potential energy of the system.";
Lagrange.help.K   = "The total kinectic energy of the system.";
if isfield(sys, 'subDiffs') && sys.subDiffs == 1

Lagrange.dL_q     = dL_q;
Lagrange.dL_dq    = dL_dq;
Lagrange.dL_dq_dt = dL_dq_dt;
Lagrange.help.dL_q     = "Cell Array of the Lagrangian differentiated w.r.t. the various GenCoords.";
Lagrange.help.dL_dq    = "Cell Array of the Lagrangian differentiated w.r.t. the various GenCoords";
Lagrange.help.dL_dq_dt = "Cell Array of the time derivatives of the Lagrangian differentiated w.r.t. the various time derivatives of the GenCoords";
end
EL = EL(t);
Lagrange.EL       = EL;
Lagrange.help.Lag = "The Lagrange equation. Lag = K-P";
Lagrange.help.EL       = "The Euler-Lagrange equations(Left hand side). EL = d/dt(dLag/d_dq) - dLad/d_q = 0";




%%% EXCTRACT MCG

if ~(isfield(sys, 'DisableExtension') && sys.DisableExtension == 1)
    
    
    
%%%%%%%%% extract Linearizable Model (M(q)*ddq + C(q,dq)*dq + G(q) = u)

[SM,EL_M,EL_C,EL_G] = ExtractMCG(EL,GenCoords) ;

Lagrange.M = EL_M;
Lagrange.C = EL_C;
Lagrange.G = EL_G;
Lagrange.help.M  = "Mass-matrix of the system model. --> M(q)*ddq + ... . MAKE SURE THIS IS INVERTIBLE!";
Lagrange.help.C  = "Coriolis-matrix of the system model. -->... + C(q,dq)*dq + ... ";
Lagrange.help.G  = "Gravity-matrix of the system model. -->... + G(q) = ...";


if isfield(sys, 'nofuncMCG') && sys.nofuncMCG == 1
    Lagrange.help.nofuncMCG = "Use this to prevent a matlabfuncton from being made(one the provides the MCG matrices)";
else
     que = [ q{:}].';
    dque = [dq{:}].';
    mass_vec = [mass{:}].';
    matlabFunction(EL_M,EL_C,EL_G,'file','MCG','vars',{que,dque,mass_vec,par});
end

end % end "extractMCG"


Lagrange.mass_ = [Lagrange.mass{:}].';
Lagrange.q_    = [Lagrange.q{:}   ].';
Lagrange.dq_   = [Lagrange.dq{:}  ].';
Lagrange.help.mass_ = "Same as 'mass', but as a vector rather than a cell array";
Lagrange.help.q_ = "Same as 'q', but as a vector rather than a cell array";
Lagrange.help.dq_ = "Same as 'dq', but as a vector rather than a cell array";


Lagrange.nq = length(Lagrange.q);
Lagrange.help.nq = "the number of states in the system";

disp('MakeLagrange --> success')
end

