function[a,b] = findParts(s,p)
% finds what two numbers, a and b, sum to s and multiply to p(product)
k = sqrt((s^2)/4-p);
a = s/2 - k;
b = s/2 + k;
end