function[out] = MakeSymbolicGenCoords(coordinates,parameters)

% This function creates sets of symbolic variables for the generalized
% coordinated and their derivatives (first and second), including one of
% thats a function of time.
% coordinates   - cell array of character strings(which are the names of the gencoords)
% parameters    - cell array of character strings(which have the parameters names)
% out           - Struct containing all the new symbolic variables



L = length(coordinates);
s = 'syms ';
sc = ';';
D = {'','d','dd'};

% Coordinates
    for i = 1:L
        name = {[coordinates{i}], [D{2} coordinates{i}], [D{3} coordinates{i}]};
        for n = 1:3
        evalin('base',[s name{n}]);
        eval([s name{n}]);
        eval(['out.' D{n} 'q{i} = ' name{n} sc]);
        end
        name = [coordinates{i} 't'];
        eval([s name '(t)']);
        eval(['out.qt{i} = ' name sc]);
    end

    
% Parameters
L = length(parameters);

    for i = 1:L
        name = {[parameters{i}]};
        evalin('base',[s name{:}]);
        eval([s name{:}]);
        eval(['out.par(i,1) = ' name{:} sc]);
    end
    
out.GenCoords = coordinates;
out.Parameres = parameters;

out.help.GenCoords = "Cells of Character-strings of variable names";
out.help.q   = "Generalized Coordinates";    
out.help.dt  = "q(t) - GenCoords as functions of time. (used in calculation of Lagrange mechanics)";
out.help.dq  = "dq/dt - Time derivatives of GenCoords";
out.help.ddq = "d2q/dt2 - Double time derivatives of GenCoords";  
out.help.par = "symbolic parameters";

evalin('base','syms g'); %gravity is always needed...

end
