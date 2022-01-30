function[A] = LeonovMatrix(f,x)

M1 = simplify(jacobian(f,x));
M2 = f*f.'/(norm(f,2));

A = simplify(M1 - M2*(M1 + M1.'));

end

