%short example of file to produce videos in matlab. 
%this video shows random noise 

% set params
fps = 30; %number of frames per seconds 
n_samples = 5 * fps; %total number of frames    

fig = figure;  %open a new figure 
mov = VideoWriter('rand.avi');  %set-up the video file 
mov.FrameRate = fps; %set the rate of frames 
open(mov);  %open the video file 

t0=0;
tf=200;
N=500;
t=linspace(t0,tf,N);    
x = sin(t/4);
y = cos(t/4);

% generate frames
for k = 1 : n_samples
%    imagesc(rand(50)); 
    xv = x(1:k);
    yv = y(1:k);
    tv = t(1:k);
    
    plot3(xv,yv,tv,'o-'); xlim([min(x),max(x)]); ylim([min(y),max(y)]); zlim([0,100])  
    %imagesc(rand(50), [0, 1]); %intensity plot of a real random matrix of size 50x50    
    drawnow; %updates the figure      
    frame = getframe(fig); %convert the figure into a frame for the video 
    writeVideo(mov,frame); %add frame to video 
end 
close(mov); %close the video file 
close(fig); %close the figure




