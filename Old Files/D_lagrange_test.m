syms q w a L m

c = cos(q);
s = sin(q);

h01 = [0 -a 0 0
       a c 0 0
       0 0 1 0
       0 0 0 1];
h12 = [1 0 0 L
       0 1 0 0
       0 0 1 0
       0 0 0 1];
h02 = h01*h12;
H=h02;

x = H(1,4);
y = H(2,4);
z = H(3,4);

vx = diff(x,q)*w;
vy = diff(y,q)*w;
vz = diff(z,q)*w;

T = 0.5*m*(vx*vx + vy*vy + vz*vz);
V = m*g*l*s;
L = T - V;

dL_dw = diff(L,w);
dL_dq = diff(L,q);
t = diff(dL_dw,q)*w + diff(dL_dw,w)*a - dL_dq;
t = simplify(t)