function[DH] = TransformationMatrix(DH)

% T - transformation from frame i-k to i
% DH - struct containing the DH-parameters for the various joints

F = fieldnames(DH.par);
M = eye(4);
M_ = M;
DH.M = M;
DH.T_.dummy = 1;
DH.A_.dummy = 1;
DH.T.dummy = 1;
DH.A.dummy = 1;
for i = 1:length(F)
    P = getfield(DH.par,F{i});
    A_ = ToNextFrame(P(1),P(2),P(3),P(4)); %Self made
    M_ = M_*A_;
    DH.T_ = setfield(DH.T_,F{i},M); % Transformation from F{i} down to F{0}
    DH.A_ = setfield(DH.A_,F{i},A_); % Transformation taking vectors from i to i-1
    %simplified versions:
    A = simplify(A_);
    M = simplify(M*A);
    DH.T = setfield(DH.T,F{i},M); % Transformation from F{i} down to F{0}
    DH.A = setfield(DH.A,F{i},A); % Transformation taking vectors from i to i-1

end

DH.M_ = M;          %Total Transformation (simplified)
DH.M = simplify(M); %Total Transformation (non simplified)


DH.help.M = "M is the total tranformation from a point in end effector frame(n-frame) and into world frame(0-frame).";
DH.help.T = "T contains the sub-total transformations. T.k -> From k to world frame.";
DH.help.A = "A contains the transformation to the previous frame. A.k -> From i to i-1.";
DH.help.A_ = "A, but not simplified";
DH.help.T_ = "T, but not simplified";
DH.help.par = "Should contain sub-fields with the parameters of the tranformation to next frame for the step it represents. MUST BE IN CORRECT ORDER: [D,TH,A,ALPHA] !!!!";

end