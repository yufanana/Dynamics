% clc
%% no damping
clear; close all
t=0:0.01:10;

m=1;
g=0;
%g=9.81;
k=5;
b=0;
L0 = 1; % original spring length
% masses = [m,m];

Fmag = 20; % magnitude of force
w = 5; % frequency units?
F = -Fmag*sin(t.*w);

Ft = 0:0.01:10; %??????

% rr = zeros(2,2); % 2x2 matrix 
initial_y1 = 0;
initial_y2 = 0;
initial_y1dot = 0;
initial_y2dot = 0;

% 1st element is the set of first order odes  
%% how to make FORCE time varying 
[t,u]=ode45( @(t,u) rhs(t,u,m,k,b,L0,g, F,Ft), t, [initial_y1 initial_y2 initial_y1dot initial_y2dot] ); 

figure; 
subplot(2,1,1); plot(t,u(:,1),'g',t,u(:,2),'r');
xlabel('t'); ylabel('y');
legend('y1','y2');
subplot(2,1,2); plot(t,F,'b');
xlabel('t'); ylabel('F');
legend('F');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dudt=rhs(t,u,m,k,b,L0,g, F,Ft) % here we write the differential equation that we want to solve   
        
        F = interp1(Ft,F,t);

        dudt = zeros(4,1);
        dudt(1) = u(3);    
        dudt(2) = u(4);
        dudt(3) = (1/m)*(F-(k/2)*(4*u(1) - 2*u(2) - 2*L0*( u(1)/(sqrt((L0^2)+u(1)^2)) - ((u(2)-u(1))/sqrt(L0^2+(u(2)-u(1))^2)) ) ) - m*g);  %y1dotdot
        dudt(4) = (1/m)*((-k/2)*(4*u(2) - 2*u(1) -2*L0*( ((u(2)-u(1))/sqrt(L0^2 + (u(2)-u(1))^2)) + u(2)/(sqrt(L0^2 + u(2)^2)))) - m*g);
end
    
