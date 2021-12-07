

function [wd_vec,zMax]=singledrivendampedoscillator_student
%computing the dynamics of a single oscillator which is driven and damped

close all;

% paramenters 
k=16; %elastic constant of spring1 
m=1; %mass of first object     

F0=1; %amplitude of forcing 
gamma=.5; %viscous resistance  

% % test of different resonance frequencies 
% wr=sqrt(k/m); %approximated resonance frequency   
% wm=wr-1; %a bit below the resonance frequency 
% wp=wr+1; %a bit above the resonance frequency 

n=200;
wd_vec=linspace(0.1,20,n);
zMax=zeros(n,1);

for i=1:n
    wd=wd_vec(i);
    param=[m,k,F0,wd,gamma]; %parameters of evolution on resonance      
    [twr,zwr]=ode45( @(t,z) rhs(t,z,param), [0 60], [0 1], param);
    
    length=size(zwr(:,1));
    final=round(0.7*length(1));
    final_z=zwr(final:length(1),1);
    zMax(i) = max(final_z);
end

figure; 
plot(wd_vec,zMax);

% param=[m,k,F0,wp,gamma]; %parameters of evolution above resonance       
% [twp,zwp]=ode45( @(t,z) rhs(t,z,param), [0 60], [0 1], param ); 
% 
% param=[m,k,F0,wm,gamma]; %parameters of evolution below resonance        
% [twm,zwm]=ode45( @(t,z) rhs(t,z,param), [0 60], [0 1], param ); 
% 
% figure; 
% plot(twm,zwm(:,1),'b',twr,zwr(:,1),'g',twp,zwp(:,1),'r'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    function dzdt=rhs(t,z,param) % here we write the differential equation that we want to solve   

        m=param(1); 
        k=param(2);         
        F0=param(3); 
        wd=param(4); 
        gamma=param(5);

        dzdt_1 = z(2);
        dzdt_2 = F0/m*sin(wd*t)-2*gamma*z(2)-(k/m)*z(1); 

        dzdt=[dzdt_1; dzdt_2];
    end     

end 







