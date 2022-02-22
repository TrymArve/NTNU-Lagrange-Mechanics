function [EL_M,EL_C,EL_G] = MCG(in1,in2,in3,in4)
%MCG
%    [EL_M,EL_C,EL_G] = MCG(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    22-Feb-2022 13:57:43

L1 = in4(2,:);
L2 = in4(3,:);
L3 = in4(4,:);
dphi = in2(2,:);
dth = in2(1,:);
dxi = in2(3,:);
g = in4(1,:);
m1 = in3(1,:);
m2 = in3(2,:);
m3 = in3(3,:);
phi = in1(2,:);
th = in1(1,:);
xi = in1(3,:);
t2 = m2+m3;
t4 = -th;
t5 = -xi;
t3 = m1+t2;
t6 = phi+t4;
t7 = phi+t5;
t8 = t5+th;
t9 = cos(t6);
t10 = cos(t7);
t11 = cos(t8);
t12 = sin(t6);
t13 = sin(t7);
t14 = sin(t8);
t15 = L2.*L3.*m3.*t10;
t16 = L1.*L3.*m3.*t11;
t17 = L1.*L2.*t2.*t9;
EL_M = reshape([L1.^2.*t3,t17,t16,t17,L2.^2.*t2,t15,t16,t15,L3.^2.*m3],[3,3]);
if nargout > 1
    EL_C = reshape([0.0,L1.*L2.*dth.*t2.*t12,-L1.*L3.*dth.*m3.*t14,-L1.*L2.*dphi.*t2.*t12,0.0,-L2.*L3.*dphi.*m3.*t13,L1.*L3.*dxi.*m3.*t14,L2.*L3.*dxi.*m3.*t13,0.0],[3,3]);
end
if nargout > 2
    EL_G = [-L1.*g.*t3.*sin(th);-L2.*g.*t2.*sin(phi);-L3.*g.*m3.*sin(xi)];
end
