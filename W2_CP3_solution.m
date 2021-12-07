%%%clean up everything 
%close all
%clear all
% //you do not have to do it ... but sometimes cleaning helps:-) 

%initializing the time vector  
t0=0;
tf=20;
N=200;
t=linspace(t0,tf,N); %//there are other ways to write vectors. which? 
dt=t(2)-t(1); 

% computing the x and y vectors 
x=3*t.^2+10*t.*cos(t);  %//note the use of ".*" for element-wise operation (that would otherwise be matrix operation)
y=t.*exp(-t/5);

%plot x vs t and y vs t 
figure; 
subplot(2,1,1); plot(t,x); xlabel('t'); ylabel('x'); %// write clear and meaningful graphs ... put labels etc.... 
subplot(2,1,2); plot(t,y); xlabel('t'); ylabel('y'); 

%3d plot 
figure; plot3(x,y,t); %// this allows you to plot a curve in 3D 
xlabel('x'); 
ylabel('y'); 
zlabel('time'); 


%compute the first derivatives in a very simple way ... how to do better?
%read the notes if you are interested  
xd=diff(x)/dt; 
yd=diff(y)/dt; 
td=t(1:end-1); 

%compute the second derivative 
xdd=diff(xd)/dt; 
ydd=diff(yd)/dt; 
tdd=t(1:N-2);

%plot of x,y, first derivatives and second derivatives     
figure; 
subplot(3,2,1); plot(t,x); xlabel('t'); ylabel('x'); %// subplot is actually pretty great to visualize more information in one slide   
subplot(3,2,2); plot(t,y); xlabel('t'); ylabel('y'); 
subplot(3,2,3); plot(td,xd); xlabel('t'); ylabel('dx/dt'); 
subplot(3,2,4); plot(td,yd); xlabel('t'); ylabel('dy/dt'); 
subplot(3,2,5); plot(tdd,xdd); xlabel('t'); ylabel('d^2 x/dt^2'); 
subplot(3,2,6); plot(tdd,ydd); xlabel('t'); ylabel('d^2 y/dt^2');



diag([1,2,3,4],3) 

%alternative way to compute derivative 
%matrix for first derivative 
D1=(-diag(ones(N,1))+diag(ones(N-1,1),1))/dt; % how was this derived? look at the notes     
xd=D1*x'; % what does the " ' " do to the vector? it is hermitian conjugation. why did I do that?  
yd=D1*y'; 
%matrix for first derivative 
D2=(diag(ones(N-1,1),-1)-2*diag(ones(N,1))+diag(ones(N-1,1),1))/dt^2; % how was this derived? look at the notes 
xdd=D2*x'; 
ydd=D2*y'; 

%plot of x,y, alternative first derivatives and second derivatives     
figure; 
subplot(3,2,1); plot(t,x); xlabel('t'); ylabel('x'); 
subplot(3,2,2); plot(t,y); xlabel('t'); ylabel('y'); 
subplot(3,2,3); plot(t(2:N-5),xd(2:N-5)); xlabel('t'); ylabel('dx/dt'); 
subplot(3,2,4); plot(t(2:N-5),yd(2:N-5)); xlabel('t'); ylabel('dy/dt'); 
subplot(3,2,5); plot(t(2:N-5),xdd(2:N-5)); xlabel('t'); ylabel('d^2 x/dt^2'); 
subplot(3,2,6); plot(t(1:N),ydd(1:N)); xlabel('t'); ylabel('d^2 y/dt^2'); 
%note that I did not plot all the points... why?   


%%%% make video  

%short example of file to produce videos in matlab. 
%this video shows random noise 

% set params
fps = 15; %number of frames per seconds 
n_samples = N; %total number of frames    

fig = figure;  %open a new figure 
mov = VideoWriter('CP3_2019.avi');  %set-up the cideo file 
mov.FrameRate = fps; %set the rate of frames 
open(mov);  %open the video file 

% generate frames
for k = 1 : n_samples
%    imagesc(rand(50)); 
    plot(x(1:k),y(1:k),'o-'); xlim([min(x),max(x)]); ylim([min(y),max(y)]);  
    %imagesc(rand(50), [0, 1]); %intensity plot of a real random matrix of size 50x50    
    drawnow; %updates the figure      
    frame = getframe(fig); %convert the figure into a frame for the video 
    writeVideo(mov,frame); %add frame to video 
end 
close(mov); %close the video file 
close(fig); %close the figure  





