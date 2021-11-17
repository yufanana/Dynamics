% clear; close all;

%% values

L_initial = 1;
% k = [2 2 2 2 2];
m = [1 1 1 1];
b = [0.01 0.01 0.01 0.01 0.01];
initial_conditions = [0 0 0 0 0 0 0 0];
Fmag = 3;

k_vec = [2 2 2 2 2;
         2 2 10 2 2;
         10 2 2 2 10;
         2 10 2 10 2];
     
k_N = size(k_vec,1);  % number of b sets to plot

w_N = 200;          % number of w points to plot
w_start = 0.1;
w_end = 5;
w_vec = linspace(w_start,w_end,w_N);

amp_vec = zeros(w_N,k_N);

%% loop

for i = 1:k_N
    
    k = k_vec(i,:);
    
    for j = 1:w_N
        w = w_vec(j);
        T = 2*pi/w;         % period
        t = 0:0.05:T*50;

        % differential equation
        [t,u] = ode45(@(t,u)rhs(t,u,m,k,b,Fmag,w), t, initial_conditions);

        start = round(0.7*size(t));       % approx region when y stabilises
        amp_max = max(u(start(1):size(t),1));
        amp_vec(j,i) = amp_max(1);
    end
end

%% plot amplitude graph
figure;  hold on;
for i = 1:k_N
    plot(w_vec,amp_vec(:,i));
end
legend( num2str(k_vec(1,:)),num2str(k_vec(2,:)),num2str(k_vec(3,:)),num2str(k_vec(4,:)) ) ;

title('Different combinations of k');
xlabel('w'); ylabel('Amplitude'); hold off;

% % natural frequencies
% xline(sqrt((3+sqrt(5))*k(1)/2*m(1)),'r');
% xline(sqrt((3-sqrt(5))*k(1)/2*m(1)),'r');
% xline(sqrt((5+sqrt(5))*k(1)/2*m(1)),'r');
% xline(sqrt((5-sqrt(5))*k(1)/2*m(1)),'r');

%% single w

% w = sqrt((3+sqrt(5))*k(1)/2*m(1));
% T = 2*pi/w;         % period
% t = 0:0.05:T*50;
% [t,u] = ode45(@(t,u)rhs(t,u,m,k,b,Fmag,w), t, initial_conditions);

%% plot individual y

% figure; 
% subplot(4,1,1); plot(t,u(:,1), 'b');
% subplot(4,1,2); plot(t,u(:,2), 'r')
% subplot(4,1,3); plot(t,u(:,3), 'g')
% subplot(4,1,4); plot(t,u(:,4), 'k')
% xlabel('t'); ylabel('y', 'Rotation',0, 'VerticalAlignment','middle','HorizontalAlignment','right');
% legend('y1','y2', 'y3', 'y4');
% 
% % F = -Fmag*sin(w*t);
% % subplot(2,1,2); plot(t,F,'b');
% % xlabel('t'); ylabel('F', 'Rotation',0, 'VerticalAlignment','middle','HorizontalAlignment','right'); 
% % legend('F');

%% video
% 
% fps = 50; %number of frames per seconds 
% 
% fig = figure;  %open a new figure 
% mov = VideoWriter('D_4m','MPEG-4');  %set-up the video file 
% mov.FrameRate = fps; %set the rate of frames 
% open(mov);  %open the video file 
% tend = T*25;
% % generate frames
% N = tend/0.05+1;
% for k = 1 : N
%     % plot mass
%     hold on; subplot(2,1,1);
%     plot(L_initial,u(k,1),'or','MarkerFaceColor','r'); hold on;
%     plot(2*L_initial,u(k,2),'ob','MarkerFaceColor','b'); hold on;
%     plot(3*L_initial,u(k,3),'ob','MarkerFaceColor','k'); hold on;
%     plot(4*L_initial,u(k,4),'ob','MarkerFaceColor','g'); hold on;
% 
%     % plot links
%     x1 = [0 L_initial];
%     y1 = [0 u(k,1)];
%     plot(x1,y1,'k');
%     x2 = [L_initial 2*L_initial];
%     y2 = [u(k,1) u(k,2)];
%     plot(x2,y2,'k'); 
%     x3 = [2*L_initial 3*L_initial];
%     y3 = [u(k,2) u(k,3)];
%     plot(x3,y3,'k'); 
%     x4 = [3*L_initial 4*L_initial];
%     y4 = [u(k,3) u(k,4)];
%     plot(x4,y4,'k'); 
%     x5 = [4*L_initial 5*L_initial];
%     y5 = [u(k,4) 0];
%     plot(x5,y5,'k');
%     hold off;
%     
%     ucat = cat(1, u(:,1), u(:,2), u(:,3), u(:,4));      % to get min/max across y1,y2,y3,y4
%     xlim([0,5*L_initial]); ylim( [min(ucat),max(ucat)] );
%     xlabel('x'); ylabel('y', 'Rotation',0, 'VerticalAlignment','middle','HorizontalAlignment','right');
%     legend('y1','y2','y3','y4');
%     
%     % plot y against time
%     subplot(2,1,2); plot(t(1:k),u(1:k,1),'r',t(1:k),u(1:k,2),'b', t(1:k),u(1:k,3),'k', t(1:k),u(1:k,4),'g');
%     xlim([0,tend]); ylim([min(ucat),max(ucat)]);
%     xlabel('t'); ylabel('y', 'Rotation',0, 'VerticalAlignment','middle','HorizontalAlignment','right');
%     legend('y1','y2','y3','y4');
% 
%     drawnow;                % updates the figure      
%     frame = getframe(fig);  % convert the figure into a frame for the video 
%     writeVideo(mov,frame);  % add frame to video 
% end 
% 
% close(mov); %close the video file 
% close(fig); %close the figure  

%% differential equation to solve 
function dudt = rhs(t,u,m,k,b,Fmag,w)
    
    y1 = u(1);
    y2 = u(2);
    y3 = u(3);
    y4 = u(4);
    y1d = u(5);
    y2d = u(6);
    y3d = u(7);
    y4d = u(8);

    m1 = m(1);
    m2 = m(2);
    m3 = m(3);
    m4 = m(4);
    
    k1 = k(1);
    k2 = k(2);
    k3 = k(3); 
    k4 = k(4);
    k5 = k(5);
    
    b1 = b(1);
    b2 = b(2);
    b3 = b(3);
    b4 = b(4);
    b5 = b(5);
    
    F = -Fmag*sin(w*t);
    
    dudt = zeros(8,1);
    dudt(1) = u(5);
    dudt(2) = u(6);
    dudt(3) = u(7);
    dudt(4) = u(8);
    dudt(5) = (1/m1)*(-k1*y1-k2*(y1-y2)-b1*(y1d)-b2*(y1d-y2d)+F);  
    dudt(6) = (1/m2)*(-k2*(y2-y1)-k3*(y2-y3)-b2*(y2d-y1d)-b3*(y2d-y3d));
    dudt(7) = (1/m3)*(-k3*(y3-y2)-k4*(y3-y4)-b3*(y3d-y2d)-b4*(y3d-y2d));
    dudt(8) = (1/m4)*(-k4*(y4-y3)-k5*(y4)-b4*(y4d-y3d)-b5*(y4d));
 
end