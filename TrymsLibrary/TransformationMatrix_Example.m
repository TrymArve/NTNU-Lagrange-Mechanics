% TransformationMatrix Example
clc;clear;close all;

% THIS IS WHAT MUST BE DONE TO USE THE FUNCTION:
syms a b c e g d1 d2 d3 th1 th2 th3
par = [a b c e g];
var = [d1 d2 d3 th1 th2 th3];

p = -(pi/2);
DH.par.A1 = [a 0 0 -p];
DH.par.A2 = [b+d1 0 0 -p];
DH.par.A3 = [c-d2 p 0 -p];
DH.par.A4 = [e+d3 0 0 0];
DH.par.A5 = [0 -th1 0 p];
DH.par.A6 = [0 th2 0 0];
DH.par.A7 = [g th3 0 0];
%%% END

DH = TransformationMatrix(DH);
DH.help


% RESULTS:

M = DH.M;

if 1
disp('A1')
disp(DH.A.A1)
disp('T1')
disp(DH.T.A1)

disp('A2')
disp(DH.A.A2)
disp('T2')
disp(DH.T.A2)

disp('A3')
disp(DH.A.A3)
disp('T3')
disp(DH.T.A3)

disp('A4')
disp(DH.A.A4)
disp('T4')
disp(DH.T.A4)

disp('A5')
disp(DH.A.A5)
disp('T5')
disp(DH.T.A5)

disp('A6')
disp(DH.A.A6)
disp('T6')
disp(DH.T.A6)

disp('A7')
disp(DH.A.A7)
disp('T7')
disp(DH.T.A7)
end

disp('M:')
disp(M)

% MAKE NUMERIC EXAMPLE:

matlabFunction(M,'file','Total_M','vars',{par,var});

P = [10 2 13 5 6];
V = [6 9 8 72 23 6];

M_example = Total_M(P,V);
disp('Example of M:')
disp(M_example)
