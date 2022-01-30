function[A] = ToNextFrame(d,th,a,alpha)
% A - The Homogeneous transformation from frame i to i+1 following the
%     DH-convention

% th    - the rotation about the z axis
% d     - the translation along z axis
% alpha - the rotation about the x axis
% a     - the translation along x axis

CLASS = string(class(th));
if CLASS == "double"
    Sth = sin(th);
    Cth = cos(th);
    if round(Cth,15) == 0
        Cth = 0;
    end
    if round(Sth,15) == 0
        Sth = 0;
    end
    
else
    Sth = sin(th);
    Cth = cos(th);
end

CLASS = string(class(alpha));
if CLASS == "double"
    Salpha = sin(alpha);
    Calpha = cos(alpha);
    if round(Calpha,15) == 0
        Calpha = 0;
    end
    if round(Salpha,15) == 0
        Salpha = 0;
    end
    
else
    Salpha = sin(alpha);
    Calpha = cos(alpha);
end

Z = [Cth -Sth 0 0;
     Sth  Cth 0 0;
      0    0  1 d;
      0    0  0 1];

X = [1    0       0   a;
     0 Calpha -Salpha 0;
     0 Salpha  Calpha 0;
     0    0      0    1];

A = Z*X;

end