function[sys] = LinearizeEL(sys,q0,dq0,masses,parameters,B)
% sys  - Struct from "MakeLagrange"
% q0   - state to linearize about(optional)
% dq0  - state to linearize about(optional)

% Form:         LM * dz       =         LCG          *  z    + LQ * U
%  -->   [I   0  ; * [dq;     = [   0          I   ; * [q ;  + [0;* U 
%         0 M(q0)]    ddq]       -dG/dq  -C(q0,dq0)]    dq]     B]

% And: dz = LA*z + LB*U  = (LM^1*LCG)*z + (LM^1*LQ)*U

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

insertpoint = 0;
if nargin > 1
    insertpoint = 1;
end

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

if insertpoint == 1
    mat = {'M', 'C', 'Gq'};
    d = {'', 'd'};
    GenCoords = sys.GenCoords; %List of variable names in q (cell array of character strings)
    Parameters = sys.par;
    Masses = sys.mass_;
    qdq = [q0(:) dq0(:)]; % concatinate q and dq


    for ii = 1:length(mat)
        
        %Substitute Variables:
        for iii = 1:length(d)
            for iiii = 1:length(GenCoords)
                text = ['syms ' d{iii} GenCoords{iiii}];
                eval(text);
    text = [mat{ii} ' = subs(' mat{ii} ',' d{iii} GenCoords{iiii} ',' num2str(qdq(iiii,iii)) ');'];
                eval(text);
            end
        end
        
        %Substitute Masses:
        for iii = 1:length(Masses)
            text = ['syms ' char(Masses(iii))];
            eval(text);
            text = [mat{ii} ' = subs(' mat{ii} ',' char(Masses(iii)) ',' num2str(masses(iii)) ');'];
            eval(text);
        end
        
        %Substitute Parameters:
        for iii = 1:length(Parameters)
            text = ['syms ' char(Parameters(iii))];
            eval(text);
            text = [mat{ii} ' = subs(' mat{ii} ',' char(Parameters(iii)) ',' num2str(parameters(iii)) ');'];
            eval(text);
        end
    end

    Gq = simplify(Gq);
    C  = simplify(C);
    M  = simplify(M);

    sys.M0  = M;  % store to struct
    sys.C0  = C;
    sys.Gq0 = Gq; %type: double (cannot simplify since its numeric)

    sys.help.M0  = "Mass matrix(M) at point of linearization";
    sys.help.C0  = "Coriolis matrix at point of linearization";
    sys.help.Gq0 = "jacobian of gravity matrix at point of linearization";

end


%Make Linear system:  ( LM*[dq;ddq] = LCG*[q;dq] + LQ*Q )
[n,m] = size(M); % DoF

LM  = blkdiag(eye(n),M); %large mass matrix
try 
    sys.invLM = LM^1;
    sys.help.invLM = "The inverse of the large Mass matrix(LM in linearization). This should always exist!";
catch
    error('OOPS!!! Large Mass matrix(LM) is not invertible!')
end
LCG = [zeros(n) eye(n);  %Large C and G matrix
         -Gq     -C  ];

if ~isfield(sys,'B')
    sys.B = eye(sys.nq); %make identity matrix of same dimention as number of coordinates(q)    
end
LQ  = [zeros(size(sys.B)); % used LQ*Q
        sys.B ];

sys.LM  = simplify(LM);     % store to struct
sys.LCG = simplify(LCG);
sys.LQ  = LQ;               % type: double (cannot simplify since its numeric)
sys.Gq  = simplify(Gq);

sys.help.LM = "Large Mass matrix - mass matrix on the left side of the linearized state space system(linearized lagrange) --> LM*[dq;ddq] = LCG*[q;dq] + LQ*Q";
sys.help.LCG = "Large C and G matrix - The matrix multiplied by [q;dq] on the right side of the linearized state space system(linearized lagrange) --> LM*[dq;ddq] = LCG*[q;dq] + LQ*Q";
sys.help.LQ = "Lagre Q - used to convert forces to the 2n state space format. [0;Q] = LQ*Q";
sys.help.Gq = "Jacobian of gravity-vector G w.r.t. q, input the linearization point q0.";

LA = simplify( LM\LCG );
LB = simplify( LM\LQ );

if insertpoint ~= 0
    disp('here')
    sys.LA = double(LA);    
    sys.LB = double(LB);
else
    sys.LA = LA;    
    sys.LB = LB;
    que = [ sys.q{:}].';
    dque = [sys.dq{:}].';
    mass = [sys.mass{:}].';
    par = sys.par;
    matlabFunction(LA,'file','LinA','vars',{que,dque,mass,par});
    matlabFunction(LB,'file','LinB','vars',{que,dque,mass,par});
end

sys.help.LA = "A matrix of Linearized EL system. ([dq;ddq] = LA*[q;dq] + LB*Q)";
sys.help.LB = "B matrix of Linearized EL system. ([dq;ddq] = LA*[q;dq] + LB*Q)";

%DISPLAY:
disp('LA:')
disp(sys.LA)
disp('LB:')
disp(sys.LB)

%CONTROLLABILITY:
if insertpoint == 1
    [Ln,Lm] = size(LA);
    Ctrb = [];
    for i = 0:Ln-1
        Ctrb = simplify( [Ctrb (LA^i)*LB] );
    end
    C_rank = rank(Ctrb);
    controllability = ((C_rank - Ln) == 0);

    sys.LA = LA;
    sys.LB = LB;
    sys.Ctrb = simplify(Ctrb);
    sys.C_rank = C_rank;
    sys.controllability = controllability;


    sys.help.Ctrb = "Controllability matrix of the Linearized state space";
    sys.help.C_rank = "Rank of the controllablility matrix";
    sys.help.controllability = "Logical: Controllability of the Linearized state space model. 1 => controllable";
    
    disp('Controllability Matrix:')
    disp(sys.Ctrb)
    disp('Rank:')
    disp(sys.C_rank)
    fprintf('Controllability: ')
    if sys.controllability
        disp('Controllable.')
    else
        disp('Not Controllable...')
    end
end




end



