function[S] = mergestructs(x,y)

% merge structs.
% use two arguments two concatinate TWO structs:
%     "mergestructs(x,y)"

% use ONE cell argument to concatinate several structs.
%     "mergestructs({x,y,z,s,t,r})"

% when several structs have equal fieldnames, the field of the LEFTmost of
% those structs will be preserved, and the others are not included in the output
% struct.


    if nargin == 1
        L = length(x);
        S = x{1};
        for ii = 2:L
            S = mergestructsTWO(S,x{ii});
        end

    else
        S = mergestructsTWO(x,y);
    end

end


function[S] = mergestructsTWO(x,y)
% Merge two structs.
% All fields of x stay untouched
% All fields of y that have DIFFERENT names than all fields of x will be
% added to x.
% (all fields of y that have the same name as a field of x, will be lost)

Cx = struct2cell(x);
Cy = struct2cell(y);

Fx = fieldnames(x);
Fy = fieldnames(y);


Y_names = string(Fy);
del = [];
for i = 1:length(Fx)
    Name = string(Fx{i});
    ind = find(Y_names == Name);
    del = [del ind];
end
Yind = 1:length(Fy);
Yind(del) = [];


S = cell2struct([Cx;Cy(Yind)],[Fx;Fy(Yind)]);

end
