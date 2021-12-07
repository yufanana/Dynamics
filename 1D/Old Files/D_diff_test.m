clear; close all;

syms y1 y2 y1d y2d L0 k1 k2 k3 m1 m2 F

d1 = sqrt(L0^2+y1^2)-L0;
d2 = sqrt(L0^2+(y2-y1)^2)-L0;
d3 = sqrt(L0^2+y2^2)-L0;

T = 0.5*(m1*y1d^2+m2*y2d^2);
V = 0.5*(k1*d1^2 + k2*d2^2 + k3*d3^2);
L = T - V;

y1dd = (1/m1)*(diff(L,y1)+F);
pretty(y1dd)
y2dd = (1/m2)*(diff(L,y2)+F);
pretty(y2dd)

% y1 = 0.5;
% y2 = 0.5;
% y1d = 0.1;
% y2d = 0.1;
% m1 = 1;
% m2 = 1;
% k1 = 1;
% k2 = 1;
% k3 = 1;
% F = 10;
% L0 = 1.5;
% 
% a1 = subs(y1dd)
% vpa(a1)