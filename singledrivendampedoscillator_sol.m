function [wd_vec,zMax]=singledrivendampedoscillator_sol
%computing the dynamics of a single oscillator which is driven and damped 

% paramenters 
k=16; %elastic constant of spring1 
m=1; %mass of first object     

F0=1; %amplitude of forcing 
gamma=.5; %viscous resistance  


% test of different resonance frequencies 
wr=sqrt(k/m); %approximated resonance frequency   
wm=wr-1; %a bit below the resonance frequency 
wp=wr+1; %a bit above the resonance frequency 

param=[m,k,F0,wr,gamma]; %parameters of evolution      
[twr,zwr]=ode45( @(t,z) rhs(t,z,param), [0 60], [0 1], param); 
param=[m,k,F0,wp,gamma]; %parameters of evolution      
[twp,zwp]=ode45( @(t,z) rhs(t,z,param), [0 60], [0 1], param ); 
param=[m,k,F0,wm,gamma]; %parameters of evolution      
[twm,zwm]=ode45( @(t,z) rhs(t,z,param), [0 60], [0 1], param ); 


figure; 
plot(twm,zwm(:,1),'b',twr,zwr(:,1),'g',twp,zwp(:,1),'r'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNCOMMENT THE FOLLOWING TO SCAN OVER THE FREQUENCIES 

% scan of many frequencies to find the resonant one 
wd_vec=linspace(.2,7,71); %frequency of forcing 
zMax=zeros(max(size(wd_vec)),1); %initialization of the vector in which we will insert the maximum amplitude of oscillations when reaching periodic regime 


for i=1:length(wd_vec) 
    
    wd=wd_vec(i); 
    %initial conditions 
    initial_x1=0;  
    initial_x1dot=1; 

    % here we do the evolution using ode45 
    % first we prepare a vector param where we write our parameters    
    param=[m,k,F0,wd,gamma];

    [t,z]=ode45( @(t,z) rhs(t,z,param), [0 500], [initial_x1 initial_x1dot], param ); 
    
    ll=max(size(t)); 
    z_eval=z(round(ll*.8):ll,1);        
    zMax(i)=max(abs(z_eval)); 
    
end;     
    
figure; plot(wd_vec,zMax,'o-'); 

  

    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    function dzdt=rhs(t,z,param) % here we write the differential equation that we want to solve   

        m=param(1); 
        k=param(2);         
        F0=param(3); 
        wd=param(4); 
        gamma=param(5); 

        dzdt_1 = z(2);    
        dzdt_2 = -k/m*z(1) + F0*cos(wd*t) - gamma*z(2);     

        dzdt=[dzdt_1; dzdt_2];
    end     

end 







