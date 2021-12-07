function W3_CP4   

% we need to solve the system of equations  
% d^2x_i/dt^2 = sum_j (x_j - x_i)*k/m_i; 
% we do a rescaling such that we can consider k/m=1 in the following (see lesson)  

close all; clear; 

t=0:0.001:10;   % time scale
% here we can simply set m1=1, m2=2 and m3=3 because we did the rescaling 
m1=1; %put some value you want to consider         
m2=2; %put some value you want to consider      
m3=3; %put some value you want to consider      

rr=rand(2,3); %generators a 2x3 matrix of random numbers    
initial_x1 = 20*(rr(1,1)-1/2); % initial x1 (randomly chosen between -10 and 10) 
initial_x2 = 20*(rr(1,2)-1/2); % initial x2 (randomly chosen between -10 and 10)  
initial_x3 = 20*(rr(1,3)-1/2); % initial x3 (randomly chosen between -10 and 10)    
initial_x1dot = 2*(rr(2,1)-1/2); % initial x1 (randomly chosen between -1 and 1)  
initial_x2dot = 2*(rr(2,2)-1/2); % initial x1 (randomly chosen between -1 and 1)  
initial_x3dot = 2*(rr(2,3)-1/2); % initial x1 (randomly chosen between -1 and 1)  

% parameter of all the masses 
masses=[m1,m2,m3]; 

% here we do the evolution using ode45 
[t,z]=ode45( @(t,z) rhs(t,z,masses), t, [initial_x1 initial_x2 initial_x3 initial_x1dot initial_x2dot initial_x3dot], masses ); 
% note the new syntax necessary to pass parameters to the function rhs 
% for example here we pass the parameters "masses"       

figure; 
plot(t,z(:,1),'-',t,z(:,2),'--',t,z(:,3),':'); % here we plot positions (indices 1, 2, 3)    
xlabel('t'); ylabel('x_1, x_2, x_3');

% total momentum: mass * velocity
totmom = m1*z(:,4) + m2*z(:,5) + m3*z(:,6); % remains constant

% total kinetic energy: sum of P^2/2M for each particle  
totkin = 0.5*m1*z(:,4).^2 + 0.5*m2*z(:,5).^2 + 0.5*m3*z(:,6).^2;

% velocity of the center of mass: momentum / mass
velcom = totmom/(m1+m2+m3);  
kincom = 0.5*(m1+m2+m3)*velcom.^2;  % kinetic energy of the center of mass 

% relative kinetic energy: v - velcom
kinrel = 0.5*m1*(z(:,4)-velcom).^2 + 0.5*m2*(z(:,5)-velcom).^2 + 0.5*m3*(z(:,6)-velcom).^2;  

% do  plots  
figure; 
subplot(2,1,1); plot(t,totmom); % here we plot the total momentum     
xlabel('t'); ylabel('Momentum of CoM');
subplot(2,1,2); plot(t, totkin,'-', t, kincom + kinrel); % here we plot the total kinetic energy and the sum of the kinetic energy of the center of mass with the relative one    
xlabel('t'); ylabel('Total K.E.');

% differential equation to solve
    function dzdt=rhs(t,z,masses)  
        dzdt_1 = z(4);    
        dzdt_2 = z(5);
        dzdt_3 = z(6);
        dzdt_4 = (z(2)+z(3)-2*z(1))/masses(1);    
        dzdt_5 = (z(3)+z(1)-2*z(2))/masses(2);
        dzdt_6 = (z(2)+z(1)-2*z(3))/masses(3);
        
        dzdt=[dzdt_1; dzdt_2; dzdt_3; dzdt_4; dzdt_5; dzdt_6];
    end

%%%% prepare the video 
% fig=figure; % make a figure 
% % now we find the minimum and maximum values of x_i so to get a good frame
% % for the figure 
% nx=[z(:,1);z(:,2);z(:,3)]; 
% mx=min(nx); 
% Mx=max(nx); 

% %%%%%%%%% now we make the video ... this is already fully done for you, but
% %%%%%%%%% you can comment it to speed-up your tests   
% % open the avi file 
% fps = 15; %number of frames per seconds 
% mov = VideoWriter('CP4_W3.avi');  %set-up the cideo file 
% mov.FrameRate = fps; %set the rate of frames 
% open(mov);  %open the video file 
% 
% N=length(t); 
% for it=1:50:N %making the movie with only one frame out of 30 
%           hold on; plot(z(it,1),t(it),'*',z(it,2),t(it),'o',z(it,3),t(it),'s'); % plotting the 3 positions individually 
%           xlim([mx Mx]); ylim([0 10]); % here we fix the frame of the figure 
%           xlabel('x_1, x_2, x_3'); ylabel('t'); % labels for the axis    
%           title('Case Problem Week 3 Prob 4');
%           frame=getframe(fig); 
%           writeVideo(mov,frame); % add the frame to the movie 
% end
% 
% close(fig);
% close(mov);



end
