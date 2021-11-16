
clear; close all

k = [2 2 2 2];
b = [1 1 1 1];
m = [1 1 1];
initial_conditions = [0 0 0 0 0 0];
Fmag = 3;

N = 200;
w_vec = linspace(0.1,5,N);
amp = zeros(N);

for j = 1:N
    w = w_vec(j);
    T = 2*pi/w;         % period
    t = 0:0.05:T*50;

    % 1st element is the set of first order odes  
    %% how to make FORCE time varying 
    [t,u]=ode45( @(t,u) rhs(t,u,m,k,b,Fmag,w), t, initial_conditions );
    
    start = round(0.7*size(t));
    amplitude = max(u(start(1):size(t),1));
    amp(j) = amplitude(1);

%     figure; 
%     subplot(3,1,1); plot(t,u(:,1));
%     xlabel('t'); ylabel('y_1');
% 
%     subplot(3,1,2); plot(t,u(:,1));
%     xlabel('t'); ylabel('y_2'); 
% 
%     subplot(3,1,3); plot(t,u(:,1));
%     xlabel('t'); ylabel('y_3');
end

% figure; 
% plot(t,u(:,1));
% xlabel('t'); ylabel('y_1');

figure; plot(w_vec,amp);
xlabel('w'); ylabel('Amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dudt = rhs(t,u,m,k,b,Fmag,w)  
        
        F = -Fmag*sin(w*t);
    
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
        
        y1 = u(1);
        y2 = u(2);
        y3 = u(3);
        y1d = u(4);
        y2d = u(5);
        y3d = u(6);
        
        dudt = zeros(6,1);
        dudt(1) = u(4);    
        dudt(2) = u(5);
        dudt(3) = u(6);
        dudt(4) = (1/m1)*(-k1*y1 - k2*(y1-y2) - b1*y1d - b2*(y1d-y2d) + F);
        dudt(5) = (1/m2)*(-k2*(y2-y1) - k3*(y2-y3) - b2*(y2d-y1d) - b3*(y2d-y3d));
        dudt(6) = (1/m3)*(-k3*(y3-y2) - k4*(y3) - b3*(y3d-y2d) - b4*(y3d));
      
end
    
