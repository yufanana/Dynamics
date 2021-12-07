% clean up everything 
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

% compute the first derivative 
xd = diff(x)/dt;
yd = diff(y)/dt;
td = t(1:end-1); 

% compute the second derivative 
xdd = diff(xd)/dt;  
ydd = diff(yd)/dt;
tdd = td(1:end-1);   

% plot of x,y, first derivatives and second derivatives     
figure; 
subplot(3,2,1); plot(t,x); xlabel('t'); ylabel('x');
subplot(3,2,2); plot(t,y); xlabel('t'); ylabel('y'); 
subplot(3,2,3); plot(td,xd); xlabel('t'); ylabel('dx/dt'); 
subplot(3,2,4); plot(td,yd); xlabel('t'); ylabel('dy/dt'); 
subplot(3,2,5); plot(tdd,xdd); xlabel('t'); ylabel('d^2 x/dt^2'); 
subplot(3,2,6); plot(tdd,ydd); xlabel('t'); ylabel('d^2 y/dt^2');

% make video  

% set params
fps = 30; %number of frames per seconds 
fig = figure;  %open a new figure 
mov = VideoWriter('HW2_Q4_video_2.avi');  %set-up the cideo file 
mov.FrameRate = fps; %set the rate of frames 
open(mov);  %open the video file 
fig.Position = [50,50,1000,700];  %specify window position and size: [left bottom width height]

% generate frames
for k = 1 : N
    
%     plot(x(1:k),y(1:k),'o-'); xlabel('x'); ylabel('y');
%     xlim([min(x),max(x)]); ylim([min(y),max(y)]);

    subplot(3,2,1); plot(t(1:k),x(1:k)); xlabel('t'); ylabel('x');
    xlim([min(t),max(t)]); ylim([min(x),max(x)]);
    
    subplot(3,2,2); plot(t(1:k),y(1:k)); xlabel('t'); ylabel('y'); 
    xlim([min(t),max(t)]); ylim([min(y),max(y)]);
    
    if k <= size(td,2)
        subplot(3,2,3); plot(td(1:k),xd(1:k)); xlabel('t'); ylabel('dx/dt'); 
        xlim([min(t),max(t)]); ylim([min(xd),max(xd)]);

        subplot(3,2,4); plot(td(1:k),yd(1:k)); xlabel('t'); ylabel('dy/dt'); 
        xlim([min(t),max(t)]); ylim([min(yd),max(yd)]);
    end
    
    if k <= size(tdd,2)
        subplot(3,2,5); plot(tdd(1:k),xdd(1:k)); xlabel('t'); ylabel('d^2 x/dt^2'); 
        xlim([min(tdd),max(tdd)]); ylim([min(xdd),max(xdd)]);

        subplot(3,2,6); plot(tdd(1:k),ydd(1:k)); xlabel('t'); ylabel('d^2 y/dt^2');
        xlim([min(tdd),max(tdd)]); ylim([min(ydd),max(ydd)]);
    end
    
    drawnow; %updates the figure      
    frame = getframe(fig); %convert the figure into a frame for the video 
    writeVideo(mov,frame); %add frame to video 
end 

close(mov); %close the video file 
close(fig); %close the figure  
