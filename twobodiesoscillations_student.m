clear; close all;
%computing normal modes for different values of masses and spring constants

% system is composed of   
% wall1 - spring1 - mass1 - spring2 - mass2 - spring3 - wall2 
% mass2 is also driven by a force F_0 sin(wt)    

% paramenters 
k1=1; %elastic constant of spring1 
m1=3; %mass of first object     
k2=1.5; %elastic constant of spring2     
m2=1; %mass of second object 
k3=0.5; %elastic constant of spring3      

wd=10; %frequency of forcing 
F0=5; %amplitude of forcing 
gamma=1; %viscous resistance  

%initial conditions 
initial_x1=0;  
initial_x2=1;  
initial_x1dot=1; 
initial_x2dot=2; 

% time scale   
t=0:0.001:100; 

%mass matrix 
M=diag([m1,m2]); 

%matrix for forces 
K=[-k1-k2, k2; k2, -k2-k3];   

%compute the normal modes and their frequencies 

[V,D]=eig(M\K); %computes the eigenvalues and eigenvectors of M^-1 * K 
% V contains the eigenvectors as its column vectors
% D is a diagonal matrix with eigenvalues along its diagonal
% in this case, eigenvalues = w0^2

w=sqrt(-diag(D)); %converts the eigenvalues to eigenfrequencies  
% minus sign because D is negative

figure; 
plot([1,2],V(:,1),'o--',[1,2],V(:,2),'o--'); 
xlim([0.5,2.5]); ylim([-1,1]);   

% dynamics
% prepare a matrix param where we write our parameters    
param = M\K;    % solves M-1 K
param_dr = [F0 gamma;wd 0];      
param = [param param_dr];

[t,z]=ode45( @(t,z) rhs(t,z,param), t, [initial_x1 initial_x2 initial_x1dot initial_x2dot], param ); 

figure;
plot(t,z(:,1),t,z(:,2)); 

% figure; 
% for i=1:500:size(t)*0.3 
%     
%     %hold on; 
%     plot([1,2],z(i,1:2),'o-','MarkerSize',20,'MarkerFaceColor','b'); 
%     xlim([0.5,2.5]); 
%     yran=max(max(abs(z(1,:))),max(abs(z(2,:))));
%     ylim([-yran,yran]); 
%     pause(0.1); 
%     
% end;     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    function dzdt=rhs(t,z,param) % here we write the differential equation that we want to solve   

        a=param(1,1); 
        b=param(1,2); 
        c=param(2,1); 
        d=param(2,2);
        % [a b F gamma; c d wdr 0]
        
        F=param(1,3); 
        wdr=param(2,3); 
        gam=param(1,4);

        dzdt_1 = z(3);    
        dzdt_2 = z(4);
        dzdt_3 = a*z(1)+b*z(2);         
        dzdt_4 = c*z(1)+d*z(2) + F*sin(wdr*t) - gam*z(4);  %here the forcing and the dissipations are added            

        dzdt=[dzdt_1; dzdt_2; dzdt_3; dzdt_4];
    end     