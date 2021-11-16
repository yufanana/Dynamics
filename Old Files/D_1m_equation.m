% inconsistent plots :/
%% no damping
% clc
close all
t=0:0.1:20;

m=1;
g=0;
k=1;
b=1;
L0 = 1; % original spring length
% masses = [m,m];

Fmag = 10; % magnitude of force
w = 4; % frequency units?
F = -Fmag*sin(t*w);
Ft = t;

% rr = zeros(2,2); % 2x2 matrix 
initial_y1 = 0;
initial_y1dot = 0;

figure; 
subplot(2,1,1); plot(t,u(:,1),'g');
xlabel('t'); ylabel('y');
legend('y');
subplot(2,1,2); plot(t,F,'b');
xlabel('t'); ylabel('F');
legend('F');


% here we do the evolution using ode45 
% 1st element is the set of first order odes  
% how to make FORCE time varying
[t,u]=ode45( @(t,u) rhs(t,u,m,k,b,L0,g, F,Ft), t, [initial_y1 initial_y1dot] ); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dudt=rhs(t,u ,m,k,b,L0,g, F,Ft) % here we write the differential equation that we want to solve   
        F = interp1(Ft,F,t); %?? why??
        
        dudt = zeros(2,1);
        dudt(1) = u(2);    
        dudt(2) = (1/m)*( F - m*g-2*k*u(1)+ (2*k*u(1)*L0/(sqrt(u(1)^2+L0^2)) ) );
end
    
