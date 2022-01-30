% Make rotation matrices(Euler angles)
function[R1,R2,R3] = MakeEulerRotationsHandle(unit)

% unit - Unit to represent angle(radians/degrees/gradians)
%        either: -number of units in 1 rotation(f.ex: pi/360/400 ...)
%            or: -one of 'deg'/'rad'/'grad' or "deg"/"rad"/"grad"


if nargin > 0
    if isnumeric(unit)
        conv = unit;
        text = ['Unit set to: ', num2str(conv),' units/rotation'];
        disp(text)
    elseif (unit == "rad") || (unit == "pi") || (unit == "radians")
        conv = pi;
        disp('Unit set to radians')
    elseif (unit == "deg") || (unit == "degrees")
        conv = 360;
        disp('Unit set to degrees')
    elseif (unit == "grad") || (unit == "gradians") || (unit == "gon")
        conv = 400;
        disp('Unit set to gradians')
    else
        conv = pi;
        disp('WARNING: Improper unit, unit set to radians.')
    end
else
    conv = pi;
end


R1 = @(phi) [1          0                 0   ;
             0 cos(phi*(pi/conv)) -sin(phi*(pi/conv));
             0 sin(phi*(pi/conv))  cos(phi*(pi/conv))];
        
R2 = @(th) [cos(th*(pi/conv))  0  sin(th*(pi/conv));
                     0         1          0        ;
           -sin(th*(pi/conv))  0  cos(th*(pi/conv))];
        
R3 = @(psi) [cos(psi*(pi/conv)) -sin(psi*(pi/conv)) 0;
             sin(psi*(pi/conv))  cos(psi*(pi/conv)) 0;
                      0                   0         1];
end

