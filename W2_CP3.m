%%%clean up everything 
close all
clear all
% //you do not have to do it ... but sometimes cleaning helps:-) 

%initializing the time vector  
t0=0;
tf=20;
N=200;
t=linspace(t0,tf,N); %//there are other ways to write vectors. which? 
dt=(tf-t0)/N;

% computing the x and y vectors 
x=3*t.^2+10*t.*cos(t);  %//note the use of ".*"  why do we use the dot? what does it do? 
y=t.*exp(-t./5); 

%plot x vs t and y vs t 
figure; 
subplot(2,1,1); plot(t,x); xlabel('t'); ylabel('x'); %// write clear and meaningful graphs ... put labels etc.... 
subplot(2,1,2); plot(t,y); xlabel('t'); ylabel('y'); 

%3d plot 
figure; plot3(x,y,t); %// this allows you to plot a curve in 3D 
xlabel('x'); ylabel('y'); zlabel('t');

%compute the first derivatives in a very simple way ... how to do better?
%read the notes if you are interested  
xd=diff(x)/dt;
yd=diff(y)/dt;
td=t(1:end-1); 

%compute the second derivative 
xdd=diff(xd)/dt;  
ydd=diff(yd)/dt;
tdd=td(1:end-1);   

%plot of x,y, first derivatives and second derivatives     
figure; 
subplot(3,2,1); plot(t,x); xlabel('t'); ylabel('x'); %// subplot is actually pretty great to visualize more information in one slide   
subplot(3,2,2); plot(t,y); xlabel('t'); ylabel('y'); 
subplot(3,2,3); plot(td,xd); xlabel('t'); ylabel('dx/dt'); 
subplot(3,2,4); plot(td,yd); xlabel('t'); ylabel('d/dty'); 
subplot(3,2,5); plot(tdd,xdd); xlabel('t'); ylabel('d^2 x/dt^2'); 
subplot(3,2,6); plot(tdd,ydd); xlabel('t'); ylabel('d^2 y/dt^2');

%alternative way to compute derivative (using self-defined matrix)
%matrix for first derivative 
D1=(-diag(ones(N,1))+diag(ones(N-1,1),1))/dt; % how was this derived? look at the notes     
XD=D1*x'; % what does the " ' " do to the vector? it is hermitian conjugation. why did I do that?  
YD=D1*y';  
%matrix for second derivative 
D2=(diag(ones(N-1,1),-1)-2*diag(ones(N,1))+diag(ones(N-1,1),1))/dt^2; % how was this derived? look at the notes 
XDD=D2*x';
YDD=D2*y'; 

%plot of x,y, alternative first derivatives and second derivatives
figure; sgtitle('Method 2') ;
subplot(3,2,1); plot(t,x); xlabel('t'); ylabel('x'); 
subplot(3,2,2); plot(t,y); xlabel('t'); ylabel('y'); 
subplot(3,2,3); plot(td,XD(1:end-1)'); xlabel('t'); ylabel('dx/dt'); 
subplot(3,2,4); plot(td,YD(1:end-1)'); xlabel('t'); ylabel('d/dty'); 
subplot(3,2,5); plot(tdd,XDD(1:end-2)'); xlabel('t'); ylabel('d^2 x/dt^2'); 
subplot(3,2,6); plot(tdd,YDD(1:end-2)'); xlabel('t'); ylabel('d^2 y/dt^2');  
%note that I did not plot all the points... why?   