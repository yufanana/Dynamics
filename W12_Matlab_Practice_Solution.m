function [w0,V,wd_vec,thetaMax,phiMax]=Matlab_Practice_Solution 


% paramenters 
g=4;  
a=1;      
gamma=0.2; 
A0=1.5; 
param=g/a; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% question 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rn=rand(1,2); 

%initial conditions 
initial_theta=pi/3*rn(1)-pi/6; %initial value of theta  
initial_phi=pi/4*rn(2)-pi/8; %initial value of phi     
initial_thetadot=0; %initial value of theta_dot    
initial_phidot=0; %initial value of phi_dot  

% time scale   
t=0:0.1:300;   

%run the ode solver   
[t,z]=ode45( @(t,z) ode_system(t,z,param), t, [initial_theta initial_phi initial_thetadot initial_phidot], param ); 


figure; 
subplot(2,1,1); plot(t,z(:,1)); ylabel('\theta'); xlabel('time');   
subplot(2,1,2); plot(t,z(:,2)); ylabel('\phi'); xlabel('time'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% questions 2 and 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%matrix of the linearized problem  
% 
Matrix=g/a*[1,-1;-1,3];    

%compute the normal modes and their frequencies 
[V,D]=eig(Matrix);     
w0=sqrt(diag(D))  

%normal modes? 
v_1=V(:,1)  
v_2=V(:,2) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  questions 4 and 5 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scanning the values of w to see where the resonances are 

wd_vec=linspace(0.1,12,50);

for i=1:length(wd_vec) 
    
    wd=wd_vec(i); 

    % here we do the evolution using ode45 
    % first we prepare a vector param2 where we write our parameters    

    param2=[g/a,gamma,A0,wd];   
    [t,z]=ode45( @(t,z) linearized_driven_ode_system(t,z,param2), t, [initial_theta initial_phi initial_thetadot initial_phidot], param2 ); 
    
    ll=max(size(t)); 
    z_eval=z(round(0.8*ll):ll,:);          
    thetaMax(i)=max(z_eval(:,1));          
    phiMax(i)=max(z_eval(:,2));   
    
end;     
    
figure; plot(wd_vec,thetaMax,'o-',wd_vec,phiMax,'x-'); 

    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    function dzdt=ode_system(t,z,param) % here we write the differential equation that we want to solve   

        goa=param; 

        dz1dt=z(3); 
        dz2dt=z(4); 
        dz3dt=-goa*sin(z(1)-z(2)); 
        dz4dt=-goa/3*((1+sin(-z(1)+3*z(2)))^3-1); 
 
        dzdt=[dz1dt,dz2dt,dz3dt,dz4dt]'; 
   

    end     



    function dzdt=linearized_driven_ode_system(t,z,param2) % here we write the linearized differential equation of the driven system that we want to solve   

        goa=param2(1);         
        gamma=param2(2);     
        A0=param2(3);       
        wd=param2(4);
        
        dz1dt=z(3); 
        dz2dt=z(4); 
        dz3dt=-gamma*z(3)-goa*(z(1)-z(2))+A0*cos(wd*t/2); 
        dz4dt=-gamma/5*z(4)-goa*(-z(1)+3*z(2)); 

        dzdt=[dz1dt,dz2dt,dz3dt,dz4dt]'; 
        
    end     


end 





