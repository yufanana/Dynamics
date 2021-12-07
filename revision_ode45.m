close all; clear;

% generate a 1x2 matrix of random numbers in the interval (0,1)
rn = rand(1,2);

% create a vector (1 row, many columns)
t_start = 0.5;
t_end = 300;
t_step = 0.1;
t = t_start:t_step:t_end;          % 0.1 is the step size

t_vec = linspace(t_start,t_end,50);

% initial conditions vector
initial_theta = pi/3*rn(1)-pi/6;  %initial value of theta  
initial_phi = pi/4*rn(2)-pi/8;    %initial value of phi     
initial_thetadot = 0;             %initial value of theta_dot    
initial_phidot = 0;               %initial value of phi_dot  
initial_conditions = [initial_theta initial_phi initial_thetadot initial_phidot];   % row vector

g=4;  
a=1;      
param=g/a; 

% ode45
[t,z] = ode45( @(t,z) ode_system(t,z,param), t, initial_conditions, param); 

% normal graph plotting
figure;
plot(t,z(:,1),t,z(:,2));
xlabel("t"); ylabel("theta");

% subplot
figure; 
subplot(2,1,1); plot(t,z(:,1));
subplot(2,1,2); plot(t,z(:,2));

% eigenvalues & natural freq
Matrix = param*[1 -1; -1 3];   
[V,D] = eig(Matrix);
diagonalD = diag(D);
w0 = sqrt(diag(D));

% for loop to get resonance
n = 100;
wd_vec = linspace(0.1,16,n);
thetaMax = zeros(n,1);
phiMax = zeros(n,1);

for i=1:n
    wd = wd_vec(i);
    param2=[g/a,0.2,1.5,wd]; 
    [t,z] = ode45(@(t,z)driven_ode_system(t,z,param2),t,initial_conditions,param2);
    ll = max(size(t));
    z_eval = z(round(0.8*ll):ll,:);
    max_values = max(z_eval);
    thetaMax(i) = max_values(1);
    phiMax(i) = max_values(2);
%     thetaMax(i) = max(z_eval(:,1));
%     phiMax(i) = max(z_eval(:,2));
end

figure;
plot(wd_vec,thetaMax,wd_vec,phiMax);
legend("Theta","Phi");

% %%%%%% create video
% tmin = min(min(x1,x2));
% tmax = max(max(x1,x2));
% zmin = min(min(y1,y2));
% zmax = max(max(y1,y2));
% 
% for k = 1 : n
%     plot(t(1:k),z(1:k,1),t(1:k),z(1:k,2)); 
%     xlabel('t'); ylabel('theta,phi'); 
%     xlim([tmin,tmax]); ylim([ymin,ymax]);
%     drawnow;                %updates the figure      
%     frame = getframe(fig);  %convert the figure into a frame for the video 
%     writeVideo(mov,frame);  %add frame to video 
% end 

%%%%%%% handy functions
% create an empty array
thetaMax = zeros(50,1);  % 50x1 array

% min / max
my_max = max(z(:,1));   % returns max element, from column 1
max_rows = max(z);      % returns row vector, with max value of each column.
A = z(:,1); B = z(:,2);
C = max(A,B);            % returns an array, with the largest elements taken from A or B.

function dzdt=ode_system(t,z,param) % here we write the differential equation that we want to solve   
        dzdt_1 = z(3);
        dzdt_2 = z(4);
        dzdt_3 = -param*sin(z(1)-z(2));
        dzdt_4 = -1/3*param*((1+sin(3*z(2)-z(1)))^3-1);
        dzdt = [dzdt_1; dzdt_2; dzdt_3; dzdt_4];    % column vector
end     

function dzdt=driven_ode_system(t,z,param2) % here we write the linearized differential equation of the driven system that we want to solve   

    goa=param2(1);         
    gamma=param2(2);     
    A0=param2(3);       
    wd=param2(4);    

    dzdt_1 = z(3);      % thetadot
    dzdt_2 = z(4);      % phidot
    dzdt_3 = -goa*(z(1)-z(2))-gamma*z(3)+A0*cos(wd*t/2);
    dzdt_4 = -goa*(3*z(2)-z(1))-gamma/5*z(4);
    dzdt = [dzdt_1; dzdt_2; dzdt_3; dzdt_4];    % column vector
end     



