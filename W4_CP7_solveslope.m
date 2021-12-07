% this program, together with traj.m, angle.m and slope_dxdt computed the
% dynamics of an object which rolls down a slope 
% the function traj.m describes the shape of the slope 
% the function angle.m returns sin(alpha)cos(alpha) where alpha is the
% angle made by the slope in a point x 
% the function slope_dxdt gives the differential dx/dt to be used in the
% ode solver. Note that x is a vector with the position and the velocity of
% the object's center of mass 

%set the parameters of simulation 
Tv=[0 8]; %time interval 
x0=[1 5]; % initial position and initial velocity 
g=9.8; % gravity 
k=1; % friction coefficient 
A=2/5; % moment of inertia divided by mass and radius square 
m=1; % mass of the object 

% solving numerically the ode 
[T,X] = ode45(@(T,X)W4_CP7_slope_dxdt(T,X,g,k,A,m),Tv,x0);

% prepare the variables xx and yy for the plot 
xx=X(:,1); 
yy=W4_CP7_traj(X(:,1)); 

% plot 
figure; plot(T,X(:,1),T,yy); xlabel('t'); ylabel('x, y');

