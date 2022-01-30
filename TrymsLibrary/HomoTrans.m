function[H] = HomoTrans(a,d,i)
% Make a Homogeneous Transformation
% a  - rotation angle [rad]
% d  - displacement distance
% id - displacement and rotation axis (which axis to drag along and rotate about)

[R1,R2,R3] = MakeEulerRotationsHandle;


switch i
    case 1
        R = R1(a);
        D = [d;0;0];
    case 2
        R = R2(a);
        D = [0;d;0];
    case 3
        R = R3(a);
        D = [0;0;d];
    otherwise
        disp('WARNING: check you arguments')
end

H = [R D; 0 0 0 1];

end
