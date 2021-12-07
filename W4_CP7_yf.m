% set the parameters of simulation 
Tv=[0 8]; %time interval 
x0=[1 5]; % initial position and initial velocity 
g=9.8; % gravity 
k=1; % friction coefficient 
A=2/5; % moment of inertia divided by mr^2
m=1; % mass of the object 

[T,X] = ode45(@(T,X)slope_dxdt(T,X,g,k,A,m), Tv, x0);

% prepare the variables xx and yy for the plot 
xx = X(:,1); 
yy = traj(X(:,1)); 

% plot 
figure; 
plot(T,xx,T,yy); 
xlabel('t'); ylabel('x, y');
legend('x','y');

function dxdt = slope_dxdt(T,X,g,k,A,m)
    sinacosa = angle(X(1));
    %dxdt = zeros(2,1);
    % X(1): x, X(@): xdot
    dxdt = [X(2); (-k/m*X(2)+g*sinacosa)/(1+A)]; 
end

function sinacosa = angle(x)
    step = 0.0001;
    % tana is the gradient
    tana = -(traj(x+step)-traj(x))/step; 
    sinacosa = tana./(1+tana.^2);
end

function y=traj(x) 
    %y=10-x;
    y=5*(cos(x)+1); 
end
    