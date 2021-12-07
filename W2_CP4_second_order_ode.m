function W2_CP4_second_order_ode

% SOLVE  d2x/dt2 - dx/dt - x = log(1+ t.^2/5)
% initial conditions: x(0) = 1, x'(0)=1 

% first thing to do is to reduce this to a system of first order
% differential equations -> 

% using x_1 = dx/dt   
% such that dx_1/dt = d^2x/dt^2

% the system of equation thus becomes 
% dx_1/dt =  x_1 + x + log(1+ t.^2/5)        
% dx/dt = x_1    


t=0:0.001:3;   % time scale

initial_y    = 1;     
initial_x    = 1; % initial position 
initial_dxdt = 1; % initial velocity   

% here we do the evolution using ode45 
[t,z]=ode45( @rhs, t, [initial_y initial_x initial_dxdt] );
% note that the initial condition is a vector. 
% this is because in this program x is a matrix with many rows and 2 columns. 
% The first columns represents the position while the second the velocity  
% x = [position, velocity]  

figure; 
subplot (3,1,1); plot(t,z(:,1)); % here we plot just the y position (index 1)    
xlabel('t'); ylabel('y');
subplot (3,1,2); plot(t,z(:,2)); % here we plot just the x position (index 2)    
xlabel('t'); ylabel('x');
subplot(3,1,3); plot(t,z(:,3)); % here we plot just the x velocity (index 3)      
xlabel('t'); ylabel('dx/dt');

    function dzdt=rhs(t,z) % here we write the differential equation that we want to solve 
        % z1=y, z2=x, z3=xdot
        dzdt = (1+t)^(-2)+t^3;  % ydot
        dzdt_1 = z(3);          % xdot
        dzdt_2 = z(3)/10 + z(2) + log(1+t^2/5);     % xdotdot
%        dxdt_2 =  - x(1); %simple test for harmonic oscillations    

        dzdt=[dzdt; dzdt_1; dzdt_2];
    end
end



