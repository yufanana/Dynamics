% example of how to use an ode solver 
% the example uses ode45. Do 
% help ode45 

close all; clear all;

%%%% first order ode %%%%%% 
% SOLVE  dy/dt = somefunction1st(t)   
% initial conditions: y(0) = 1 
% the ode solver returns y(t) for certain values of t  

t=0:0.1:10;   % points of evaluation of the function y       
initial_y=1;  % initial condition y(0) = 1   

%here you call ode45 to do the numerical integration for you ... you can choose between various ode solvers  
%type help ode45 
[t,y]=ode45( @somefunction1st, t, initial_y);  

figure; 
plot(t,y,'o',t,5*exp(-t)-4*ones(101,1),'r-'); % plot the numerical results with some circles and compare them to the analytical function drawn by a red line        
xlabel('t'); %sets the y-axis label     
ylabel('y'); %sets the y-axis label   

%%%% second order ode %%%%%% 
% SOLVE  d^2 y / dt^2 = somefunction2nd(t)   
% initial conditions: y(0) = 1   and  y'(0) = 1 
% the ode solver returns y(t) for certain values of t  

t=0:0.1:10;   % points of evaluation of the function y       
initial_y=1;  % initial condition y(0) = 1  
initial_yp=0;  % initial condition y'(0) = 1  

%here you call ode45 to do the numerical integration for you ... you can choose between various ode solvers  
%type "help ode45" 
[t,y]=ode45( @somefunction2nd, t, [initial_y, initial_yp]);  

figure; 
subplot(2,1,1); plot(t,y(:,1),'o',t,cos(t),'r-'); % note that we only select one column of y because the other one contains the velocities   
xlabel('t'); %sets the y-axis label     
ylabel('y'); %sets the y-axis label   
subplot(2,1,2); plot(t,y(:,2),'o',t,-sin(t),'r-'); % note that we only select one column of y because the other one contains the velocities   
xlabel('t'); %sets the y-axis label     
ylabel('y'''); %sets the y-axis label 

