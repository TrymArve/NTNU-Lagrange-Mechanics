function[expr] = RepairInertia(expr,from,to)
% This function is used to rebuild the inertia model of a mechanical system.
% When using the "MakeLagrange" function to build your euler-lagrange
% equations, you may only have point masses. Therefore you cannot implement
% inertia directly. To circumvent this issue, you may use three points for
% every mass. One to represent the actual mass-point, and then two points in
% opposite directions from the main point. Their raduis MUST BE: R = 1 in
% each direction. Their average lies in the main point, and does therefore
% not impact the translational mechanics, but introduce inertia. 

% To rebuild, simply make the two masses each equal to half the inertia:
% m11,m12 = (1/2)*I1, such that their sum is I1. Then change the mass of
% the main point to the actual mass minus the minus the inertia(the
%    remaining mass, that is not added to the middle my the two other masses)

%EXAMPLE: 
% - A point has mass M.
% - You model it as three point masses: m1(main) , m11 , m12
% - - Positions: 
% - - - m1: pos_1
% - - - m11: pos_11 = pos_1 + [cos(q1);sin(q1)]
% - - - m12: pos_12 = pos_1 - [cos(q1);sin(q1)]
% - - - (where q1 is the rotation of the mass)
% - - sys.pos = {pos_1, pos_11, pos_12}
% - Change the masses of the resulting EL,M,C,G,W...
% - - Change: m11 --> (1/2)*I1 , m12 --> (1/2)*I1 , m1 --> (m1 - I1)

% FUNCTION:
% expr - the expression to repair
% from - cells containing the variable names to change(the old ones)
%        (always three masses: from = {'m1','m11','m12'} )
% to   - Cells containing the mass and ineria variable name to change to(new)
%        (always two variables, mass and inertia for the point: to = {'I1','m1'} )



f1 = from{1}; %main mass
f2 = from{2}; %peripheral mass 1
f3 = from{3}; %peripheral mass 2

if length(to) > 1
    m = to{1}; %mass
else
    m = f1;
end
I = to{2}; %inertia


s = ' ';
eval(['syms ' f1 s f2 s f3 s I s m]);
evalin('base',['syms ' I s m]);

before = 'expr = subs(expr,';
mid = ',(1/2)*';
after = ');';
text1 = [before f2 mid I after];
text2 = [before f3 mid I after];
eval(text1);
eval(text2);


if to{1} == 0
    text3 = [before f1 ',0' after];  %if mass is in origin --> only inertia
else
    text3 = [before f1 ',(' m '-' I ')' after];
end
eval(text3);

% expr = subs(expr,m2,(1/2)*I1);
% expr = subs(expr,m3,(1/2)*I1);
% expr = subs(expr,m1,(m1-I1));

expr = simplify(expr);

end


% %% Repair Inertia(old version)
% syms m1 m2 m3 m4 m5 m6 I1 I2
% EL = S.EL;
% EL = subs(EL,m2,(1/2)*I1);
% EL = subs(EL,m3,(1/2)*I1);
% EL = subs(EL,m5,(1/2)*I2);
% EL = subs(EL,m6,(1/2)*I2);
% EL = subs(EL,m1,(m1-I1));
% EL = subs(EL,m4,(m2-I2));
% EL = simplify(EL);
