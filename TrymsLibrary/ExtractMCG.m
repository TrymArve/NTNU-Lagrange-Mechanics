function[SM,M,C,G] = ExtractMCG(EM,q,makefunc)
% EM - equations of motion (excluding control forces Q)
% form: EM = [0; 0; ...]
% q - cells containing character arrays of the variable names

func = 0;
if nargin > 2
    func = makefunc;
end



%%%%%%%%% SET SYMBOLIC VARIABLES
    n = length(q);

    s = 'syms ';
    sc = ';';
    for i = 1:n
        text = [s 'd' q{i} sc];
        eval(text)
        text = ['dq{' num2str(i) '} = d' q{i} sc];
        eval(text)
        text = [s 'dd' q{i} sc];
        eval(text)
        text = ['ddq{' num2str(i) '} = dd' q{i} sc];
        eval(text)
    end
%%%%%%%%
 q   = [  q{:}].';
 dq  = [ dq{:}].';
 ddq = [ddq{:}].';

%Find mass matrix
M = simplify( jacobian(EM,ddq) );

 
%Find Coriolis matrix
C = [];
for i = 1:n
    H = hessian(EM(i),dq);
    T = EM(i) - (1/2)*dq.'*H*dq;
    jT = jacobian(T,dq);
    C_vec = jT + (1/2)*dq.'*H;
    C = simplify([C; C_vec ]);
end


% extract REMANING TERMS
G = simplify( EM - M*ddq - C*dq );

% Rebuild system model for comparison to EoM
SM = simplify( M*ddq + C*dq + G );

if func
    
    %Find Parameters:
    par = symvar(SM)
    %Remove gencoords:
    ind = par - q(:)
    [ind1,ind2] = find(ind == 0)
    par(ind2) = []
    %Remove gencoords derivatives:
    ind = par - dq(:)
    [ind1,ind2] = find(ind == 0)
    par(ind2) = []
    % Remove double derivatives:
    ind = par - ddq(:)
    [ind1,ind2] = find(ind == 0)
    par(ind2) = []

    matlabFunction(M,C,G,'file','MCG','Vars',{q,dq})
end

end