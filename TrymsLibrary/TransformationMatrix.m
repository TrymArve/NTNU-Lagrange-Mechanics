function[DH] = TransformationMatrix(DH)

% T - transformation from frame i-k to i
% DH - struct containing the DH-parameters for the various joints

F = fieldnames(DH.par);
M = eye(4);
DH.M = M;
DH.T.dummy = 1;
DH.A.dummy = 1;
for i = 1:length(F)
    P = getfield(DH.par,F{i});
    A = ToNextFrame(P(1),P(2),P(3),P(4)); %Self made
    M = M*A;
    DH.T = setfield(DH.T,F{i},M); % Transformation from F{i} down to F{0}
    DH.A = setfield(DH.A,F{i},A); % Transformation taking vectors from i to i-1
end

DH.M = M; %Total Transformation

end