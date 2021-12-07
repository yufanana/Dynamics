clear;
% ch: added gravity, damping

%% obtain symbolic expression for y1dd & y2dd
syms y1 y2 y3 y1d y2d y3d L0 k1 k2 k3 k4 b1 b2 b3 b4 m1 m2 m3 g F(t)

% spring extensions
d1 = sqrt(L0^2+y1^2)-L0;
d2 = sqrt(L0^2+(y2-y1)^2)-L0;
d3 = sqrt(L0^2+(y3-y2)^2)-L0;
d4 = sqrt(L0^2+y3^2)-L0;

% individual kinetic energies
T1 = 0.5*m1*y1d^2;
T2 = 0.5*m2*y2d^2;
T3 = 0.5*m3*y3d^2;
T = T1+T2+T3;

% individual spring energies
Vk1 = 0.5*k1*d1^2; 
Vk2 = 0.5*k2*d2^2;
Vk3 = 0.5*k3*d3^2;
Vk4 = 0.5*k4*d4^2;
Vk = Vk1+Vk2+Vk3+Vk4;

% individual gravitational energies
Vg1 = g*m1*y1;
Vg2 = g*m2*y2;
Vg3 = g*m3*y3;
Vg = Vg1+Vg2+Vg3;

%V = Vk+Vg;

L = T - Vk;

% fb damping force
fb1 = -b1*y1d;
fb2 = -b2*(y2d-y1d);
fb3 = -b3*(y3d-y2d);
fb4 = -b4*y3d;

% non-conservative force
Wb1 = -b1*(y1d-0)*(y1-0);
Wb2 = -b2*(y2d-y1d)*(y2-y1);
Wb3 = -b3*(y3d-y2d)*(y3-y2);
Wb4 = -b4*(0-y3d)*(0-y3);
Wb = Wb1+Wb2+Wb3+Wb4;

Q1 = F(t)+diff(Wb,y1); 
Q2 = diff(Wb,y2);
Q3 = diff(Wb,y3);

% y = zeros(3,1);
y1dd = (1/m1)*(diff(L,y1)+Q1);
y2dd = (1/m2)*(diff(L,y2)+Q2);
y3dd = (1/m3)*(diff(L,y3)+Q3);

%% values
g = 0;

L_initial = 1;
k = [2 2 2 2];
b = [0.1 0.1 0.1 0.1];
m = [1 1 1];
initial_conditions = [0 0 0 0 0 0];
Fmag = 3;

w = 0.25;
T = 2*pi/w;         % period
t = 0:0.05:T*15;

%% differential equation
[t,u] = ode45(@(t,u)rhs(t,u,m,k,b,L_initial,Fmag,w,y1dd,y2dd,y3dd), t, initial_conditions);

%% plot results
figure; 
subplot(2,1,1); plot(t,u(:,1), 'b',t,u(:,2),'r',  t,u(:,3),'k');
xlabel('t'); ylabel('y', 'Rotation',0, 'VerticalAlignment','middle','HorizontalAlignment','right');
legend('y1','y2', 'y3');
title(['w = ' num2str(w)]);

F = -Fmag*sin(w*t);
subplot(2,1,2); plot(t,F,'b');
xlabel('t'); ylabel('F', 'Rotation',0, 'VerticalAlignment','middle','HorizontalAlignment','right'); 
legend('F');

%% differential equation to solve 
function dudt = rhs(t,u,m,k,b,L_initial,Fmag,w,y1dd,y2dd,y3dd)
    
    % give values to the symbols
    F = -Fmag*sin(w*t);
    L0 = L_initial;
    
    y1 = u(1);
    y2 = u(2);
    y3 = u(3);
    y1d = u(4);
    y2d = u(5);
    y3d = u(6);
    
    m1 = m(1);
    m2 = m(2);
    m3 = m(3);
    k1 = k(1);
    k2 = k(2);
    k3 = k(3); 
    k4 = k(4);
    b1 = b(1);
    b2 = b(2);
    b3 = b(3);
    b4 = b(4);

    % evaluate the value of symbolic expression for y1dd & y2dd  
    y1dd = subs(y1dd);              % symbolic substitution
    y2dd = subs(y2dd);
    y3dd = subs(y3dd);
    y1dd = double(vpa(y1dd));       % variable-precision arithmetic, convert to double
    y2dd = double(vpa(y2dd));
    y3dd = double(vpa(y3dd));

    dudt = [y1d; y2d; y3d; y1dd; y2dd; y3dd];
end
