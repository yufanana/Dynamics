close all;

% initialize the time vector  
t0 = 0;
tf = 3;
N = 200;
t = linspace(t0,tf,N);
dt = (tf-t0)/N;

% compute the x and y vectors 
y = exp(-t.^4/6);
x = ((1+t).^2).^-1 + y.*sin(t);

%%% NOTE: diff reduces array size by 1 

% compute the first derivative
xd = diff(x)/dt;
yd = diff(y)/dt;
td = t(1:end-1); 

% compute the second derivative 
xdd = diff(xd)/dt;  
ydd = diff(yd)/dt;
tdd = td(1:end-1);   

% element-wise multiplication
posmod = 5*sqrt(x.^2 + y.^2);
velmod = 5*sqrt(xd.^2 + yd.^2);

% 2d plot
figure; 
subplot(2,1,1); plot(t,x); xlabel('t'); ylabel('x'); %// write clear and meaningful graphs ... put labels etc.... 
subplot(2,1,2); plot(t,y); xlabel('t'); ylabel('y'); 

% 3d plot 
figure; plot3(x,y,t);
xlabel('x'); ylabel('y'); zlabel('t');
