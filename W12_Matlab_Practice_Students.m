function [w0,V,wd_vec,thetaMax,phiMax]=W12_Matlab_Practice_Students       

close all;

% parameters 
g=4;  
a=1;      
gamma=0.2; 
A0=1.5; 
param=g/a; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% question 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initial conditions 
rn=rand(1,2); 
initial_theta=pi/3*rn(1)-pi/6;  %initial value of theta  
initial_phi=pi/4*rn(2)-pi/8;    %initial value of phi     
initial_thetadot=0;             %initial value of theta_dot    
initial_phidot=0;               %initial value of phi_dot  
initial_conditions = [initial_theta initial_phi initial_thetadot initial_phidot]; 

% time scale   
t=0:0.01:300;   

%run the ode solver   
% [t,z] = ode45(@(t,z) ode_system(t,z,param),t,initial_conditions, param);
[t,z]=ode45( @(t,z) ode_system(t,z,param), t, initial_conditions, param ); 

figure; 
subplot(2,1,1); plot(t,z(:,1)); xlabel("t"); ylabel("\theta");
subplot(2,1,2); plot(t,z(:,2));  xlabel("t"); ylabel("\Phi");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% questions 2 and 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%matrix of the linearized problem  

Matrix = param*[1 -1; -1 3];   

%compute the normal modes and their frequencies 
[V,D]=eig(Matrix);    

w0 = sqrt(diag(D));

%normal modes? 
%each column is the eigenvector
figure;
plot([1,2],V(:,1),'o--',[1,2],V(:,2),'o--'); 
xlim([0.5,2.5]); ylim([-1,1]);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  questions 4 and 5 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scanning the values of w to see where the resonances are 

n = 100;
wd_vec = linspace(0.1,16,n);
thetaMax = zeros(n,1);
phiMax = zeros(n,1);

for i=1:length(wd_vec) 
    
    wd=wd_vec(i); 

    % here we do the evolution using ode45 
    % first we prepare a vector param2 where we write our parameters    

    param2=[g/a,gamma,A0,wd];   
    [t,z]=ode45( @(t,z) linearized_driven_ode_system(t,z,param2),t,initial_conditions,param2 ); 
    
    ll = max(size(t)); 
    z_eval= z(round(0.8*ll):ll,:);
    thetaMax(i)= max(z_eval(:,1)); 
    phiMax(i) = max(z_eval(:,2));

end   
    
figure; 
plot(wd_vec,thetaMax,'o-',wd_vec,phiMax,'x-'); 
legend("\theta Max","\phi Max");
xlabel("wd"); ylabel("Max");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    function dzdt=ode_system(t,z,param) % here we write the differential equation that we want to solve   
        dzdt_1 = z(3);
        dzdt_2 = z(4);
        dzdt_3 = -param*sin(z(1)-z(2));
        dzdt_4 = -1/3*param*((1+sin(3*z(2)-z(1)))^3-1);
        dzdt = [dzdt_1; dzdt_2; dzdt_3; dzdt_4];    % column vector
    end     

    function dzdt=linearized_driven_ode_system(t,z,param2) % here we write the linearized differential equation of the driven system that we want to solve   
        
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


end 





