function[out] = MakeSymbolicGenCoords(in)

% This function creates sets of symbolic variables for the generalized
% coordinated and their derivatives (first and second), including one of
% thats a function of time.
% in  - cell array of character strings(which are the names of the gencoords)
% out - Struct containing all the new symbolic variables

L = length(in);
s = 'syms ';
sc = ';';
D = {'','d','dd'};

    for i = 1:L
        name = {[in{i}], [D{2} in{i}], [D{3} in{i}]};
        for n = 1:3
        evalin('base',[s name{n}]);
        eval([s name{n}]);
        eval(['out.' D{n} 'q{i} = ' name{n} sc]);
        end
        name = [in{i} 't'];
        %evalin('base',[s name '(t)']);
        eval([s name '(t)']);
        eval(['out.qt{i} = ' name sc]);
    end

out.GenCoords = in;

out.help.GenCoords = "Cells of Character-strings of variable names";
out.help.q   = "Generalized Coordinates";    
out.help.dt  = "q(t) - GenCoords as functions of time. (used in calculation of Lagrange mechanics)";
out.help.dq  = "dq/dt - Time derivatives of GenCoords";
out.help.ddq = "d2q/dt2 - Double time derivatives of GenCoords";  

evalin('base','syms g');


%disp('success!!')
end
   


%%%%%%%%%%%%%%% OLD

% % States (required)
% syms x1t(t) x2t(t) e1t(t) e2t(t) % (MUST be functions of time)
% mysystem.qt = {x1t, x2t, e1t, e2t};
% 
% 
% % Derivatives  
% %(the indecies must be consistent)
% %(i'th element dq{i} is derivative of q{i}, etc.)
% % States (required)
% syms x1 x2 e1 e2 % (Not functions of time)
% mysystem.q = {x1, x2, e1, e2};
% syms dx1 dx2 de1 de2 ddx1 ddx2 dde1 dde2 %(NOT functions of time)
% mysystem.dq =  { dx1,  dx2,  de1,  de2}; 
% mysystem.ddq = {ddx1, ddx2, dde1, dde2};