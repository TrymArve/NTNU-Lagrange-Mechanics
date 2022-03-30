function[Lin] = LinearizeModel(sys)
% Use to linearize the model:
% (M(q)*ddq + C(q,dq)*dq + G(q) = Q)
% to:
% dr = A(r0)*r + r0 + B(U=0)*U
% where: r = [q; dq], U = Q(non-zero)

M    = sys.M;
C    = sys.C;
G    = sys.G;
q    = [sys.q{:}].';
dq   = [sys.dq{:}].';
ddq  = [sys.ddq{:}].';
Q    = [sys.Q{:}].';
par  = sys.par;
mass = [sys.mass{:}].';

Lin = sys;

%Define U
U = Q(find(Q~=0));
Lin.U = U;


% find ddq dynamics: (ddq = N(q,dq,u))
N = simplify( M\(Q - C*dq - G) );
Lin.N = N;

A_nl = [dq; N];

% Linearize to get on form (dr = A(r0)*r + r0 + B(Q=0)*Q)
r = [q;dq];
A = simplify( jacobian(A_nl,r) );
B = simplify( jacobian(A_nl,U) );
Lin.A = A;
Lin.B = B;

if isfield(sys, 'nofunc') && sys.nofunc == 1
else
matlabFunction(A,B,'file','AB','vars',{r,U,mass,par});
end


end