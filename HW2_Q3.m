
function HW2_Q3

% clean up 
close all;

% time scale
t0 = 0;
tf = 1.4; % gives an error for t>1.4
N = 200;
t=linspace(t0,tf,N);
dt = (tf-t0)/N;

initial_x    = 1; 
initial_z    = 1;
initial_dxdt = 0;
initial_dzdt = 2;


% compute the second derivative 
xdd = diff(u(:,3))/dt;  
zdd = diff(u(:,4))/dt;
tdd = t(1:end-1);   

% compute modulus
posmod = sqrt(u(:,1).^2 + u(:,2).^2);
velmod = sqrt(u(:,3).^2 + u(:,4).^2);
accmod = sqrt(xdd.^2 + zdd.^2);

% display modulus results
figure; 
subplot (3,1,1); plot(t,posmod);    
xlabel('t'); ylabel('Position Mod');
subplot (3,1,2); plot(t,velmod);    
xlabel('t'); ylabel('Velocity Mod');
subplot (3,1,3); plot(tdd,accmod);    
xlabel('t'); ylabel('Acceleration Mod');

% display vector results
figure; 
subplot (6,1,1); plot(t,u(:,1)); % plot the x position (index 1)    
xlabel('t'); ylabel('x');
subplot (6,1,2); plot(t,u(:,2)); % plot the z position (index 2)    
xlabel('t'); ylabel('z');
subplot (6,1,3); plot(t,u(:,3)); % plot the x velocity (index 3)      
xlabel('t'); ylabel('dx/dt');
subplot (6,1,4); plot(t,u(:,4)); % plot the z velocity (index 4)      
xlabel('t'); ylabel('dz/dt');
subplot (6,1,5); plot(tdd,xdd); % plot the x acceleration  
xlabel('t'); ylabel('d^2x/dt^2');
subplot (6,1,6); plot(tdd,zdd); % plot the z acceleration   
xlabel('t'); ylabel('d^2z/dt^2');

% make video  

% set params
fps = 30; %number of frames per seconds 

fig = figure;  %open a new figure 
mov = VideoWriter('HW2_Q3_video.avi');  %set-up the video file 
mov.FrameRate = fps; %set the rate of frames 
open(mov);  %open the video file 

% generate frames
for k = 1 : N
    plot(u(1:k,1),u(1:k,2),'o-'); xlabel('x'); ylabel('z');
    xlim([min(u(:,1)),max(u(:,1))]); ylim([min(u(:,2)),max(u(:,2))]);

    drawnow; %updates the figure      
    frame = getframe(fig); %convert the figure into a frame for the video 
    writeVideo(mov,frame); %add frame to video 
end 

close(mov); %close the video file 
close(fig); %close the figure  

end