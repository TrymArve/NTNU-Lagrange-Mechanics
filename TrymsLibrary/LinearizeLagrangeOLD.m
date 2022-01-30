function[sys] = LinearizeLagrangeOLD(sys,q0,dq0,par,mass)
% sys  - Struct from "MakeLagrange"
% q0   - state to linearize about
% dq0  - state to linearize about
% par  - parameter values
% mass - mass values

%(remember to update sys.mass/sys.q/sys.dq after repairing inertia!!!)

% MUST HAVE:
% (sys.q    - Symbolic variables)
% (sys.dq   - Symbolic variables)
% (sys.par  - Symbolic parameters)
% (sys.mass - Symbolic masses)
% (sys.B    - Numeric control force matrix: Mddq + Cdq + G = B*U)

% AND:
% (sys.EL        - Euler-Lagrange equations(to linearize))
% (sys.GenCoords - Cells containing names of variables)
% OR:
% (sys.M  - Symbolic mass matrix)
% (sys.C  - Symbolic coriolis matrix)
% (sys.G  - Symbolic gravity matrix)

%Extract Matrices
if isfield(sys,'EL') && isfield(sys,'GenCoords')
    [EL,M,C,G] = ExtractMCG(sys.EL,sys.GenCoords);
elseif isfield(sys,'M') && isfield(sys,'C') && isfield(sys,'G')
      M = sys.M;
      C = sys.C;
      G = sys.G;
else
    error('OOPS!! Could not obtain matrix form of EL equations.')
end

%Linearize G:
Gq = jacobian(G,[sys.q{:}].'); % Linearize wrt q (not dependent on dq)

%Find Numeric Matrices:
matlabFunction(M,C,Gq,'file','MCG_Lin','vars',{[sys.q{:}].',[sys.dq{:}].',sys.par,sys.mass});
[M,C,Gq] = MCG_Lin(q0,dq0,par,mass);

%Make Linear system:  ( LM*[dq;ddq] = LCG*[q;dq] + LQ*Q )
n = length(q0); % DoF

LM  = blkdiag(eye(n),M); %large mass matrix
try 
    sys.invLM = LM^1;
catch
    error('OOPS!!! Large Mass matrix(LM) is not invertible!')
end
LCG = [zeros(n) eye(n);  %Lagra C and G matrix
         -Gq     -C  ];
LQ  = [zeros(n); % used LQ*Q
        sys.B ];

sys.LM  = LM;    % store to struct
sys.LCG = LCG;
sys.LQ  = LQ;
sys.Gq  = Gq;

sys.help.LM = "Large Mass matrix - mass matrix on the left side of the linearized state space system(linearized lagrange) --> LM*[dq;ddq] = LCG*[q;dq] + LQ*Q";
sys.help.LCG = "Large C and G matrix - The matrix multiplied by [q;dq] on the right side of the linearized state space system(linearized lagrange) --> LM*[dq;ddq] = LCG*[q;dq] + LQ*Q";
sys.help.LQ = "Lagre Q - used to convert forces to the 2n state space format. [0;Q] = LQ*Q";
sys.help.Gq = "Jacobian of gravity-vector G w.r.t. q, input the linearization point q0.";

LA = LM\LCG;
LB = LM\LQ;
Ctrb = ctrb(LA,LB);
C_rank = rank(Ctrb);
controllability = ((C_rank - 2*n) == 0);

sys.LA = LA;
sys.LB = LB;
sys.Ctrb = Ctrb;
sys.C_rank = C_rank;
sys.controllability = controllability;

sys.help.LA = "A matrix of Linearized EL system. ([dq;ddq] = LA*[q;dq] + LB*Q)";
sys.help.LB = "B matrix of Linearized EL system. ([dq;ddq] = LA*[q;dq] + LB*Q)";
sys.help.Ctrb = "Controllability matrix of the Linearized state space";
sys.help.C_rank = "Rank of the controllablility matrix";
sys.help.controllability = "Logical: Controllability of the Linearized state space model. 1 = controllable";



end



