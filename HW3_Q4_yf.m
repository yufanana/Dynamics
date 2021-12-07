close all; clear;

N = 200;
t0 = 0;
tf = 5;
t = linspace(t0,tf,N);
k = 1;
m1 = 1;
m2 = 2;

width = 5;
initial_r1 = ([rand; rand]-0.5).*width;
initial_r2 = ([rand; rand]-0.5).*width;
initial_r1dot = ([rand; rand]-0.5).*width;
initial_r2dot = ([rand; rand]-0.5).*width;

% initial_r1 = [1; 2];
% initial_r2 = [-5; 5];
% initial_r1dot = [-3; 3];
% initial_r2dot = [3; 7];

initial_conditions = [initial_r1(1) initial_r1dot(1) initial_r2(1) initial_r2dot(1) initial_r1(2) initial_r1dot(2) initial_r2(2) initial_r2dot(2)];

[t,u] = ode45( @(t,u)rhs(t,u,k,m1,m2), t, initial_conditions);

% angular momentum about origin
x1 = u(:,1);
x1dot = u(:,2);
x2 = u(:,3);
x2dot = u(:,4);
y1 = u(:,5);
y1dot = u(:,6);
y2 = u(:,7);
y2dot = u(:,8);

H1 = (x1.*y1dot - y1.*x1dot).*m1;
H2 = (x2.*y2dot - y2.*x2dot).*m2;
totalH = H1+H2; % should be constant over time

% plot
figure;
subplot(3,1,1); plot(u(:,1),u(:,5)); 
xlabel('x1'); ylabel('y1');
subplot(3,1,2); plot(u(:,3),u(:,7)); 
xlabel('x1'); ylabel('y2');
subplot(3,1,3); plot(t,totalH); 
xlabel('t'); ylabel('total H');

% make video  

% set params
fps = 15; %number of frames per seconds 
fig = figure;  %open a new figure 
mov = VideoWriter('HW3_Q4_video.avi');  %set-up the cideo file 
mov.FrameRate = fps; %set the rate of frames 
open(mov);  %open the video file 
fig.Position = [50,50,1000,700];  %specify window position and size: [left bottom width height]

xmin = min(min(x1,x2));
xmax = max(max(x1,x2));
ymin = min(min(y1,y2));
ymax = max(max(y1,y2));

% generate frames
for k = 1 : N
    
    plot(x1(1:k),y1(1:k), x2(1:k),y2(1:k)); 
    xlabel('x1,x2'); ylabel('y1,y2'); 
    xlim([xmin,xmax]); ylim([ymin,ymax]);
    
    drawnow; %updates the figure      
    frame = getframe(fig); %convert the figure into a frame for the video 
    writeVideo(mov,frame); %add frame to video 
end 

close(mov); %close the video file 
close(fig); %close the figure  


function dudt = rhs(t,u,k,m1,m2)
    r1 = [u(1); u(5)];  % r1 = [x1; y1]
    r2 = [u(3); u(7)];  % r2 = [x2; y2]
    f = forces(k,r1,r2);
    % f = [f12x, f12y, f21x, f21y];
    
    dudt_1 = u(2);
    dudt_2 = f(3)/m1;
    dudt_3 = u(4);
    dudt_4 = f(1)/m2;
    dudt_5 = u(6);
    dudt_6 = f(4)/m1;
    dudt_7 = u(8);
    dudt_8 = f(2)/m2;

    dudt = [dudt_1; dudt_2; dudt_3; dudt_4; dudt_5; dudt_6; dudt_7; dudt_8];
end

function f = forces(k,r1,r2)
    % fij is interpreted as force of i on j
    f12 = (norm(r1-r2)^2).*(r1-r2).*k;  % norm returns magnitude of the vector
    f21 = (norm(r1-r2)^2).*(r2-r1).*k;
    % f = [f12x, f12y, f21x, f21y];
    f = [f12(1), f12(2), f21(1), f21(2)];
end
