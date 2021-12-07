%%%%%%%%%
% vectors
%%%%%%%%%

a=[1,3,5,7]         % row vector
b=[1;3;5;7]         % column vector
b'                  % transpose
c=(1:2:7)           % (start,step,end)
d=linspace(1,7,4)   % linearly spaced vectors: (start,end,#points)
e=logspace(1,7,4)   % log spaced vectors: (start,end,#points) 

%%%%%%%%%%
% matrices
%%%%%%%%%%

A=[1,2;3,4]  
A' 
B=[1,2,3;4,5,6]  
C=ones(2,4)     % ones(rows,cols)
D=ones(3)  
E=zeros(2,3)  
F=zeros(2)          % square matrix
G=rand(2,4)         % unif dist (start,end)
H=diag([1,2,3,4])   % creates square matrix
HH =diag([1,2,3,4],1)   % adds diagonal on column 1  (5x5)
I=diag(H)          % gets the diag values
J=eye(3)          % identity matrix

%%%%%%%
% plots
%%%%%%%

x=linspace(-10,10,201); 
y=x.^2+3-x;     % . to denote element-wise operation, instead of matrix operation
figure 
plot(x,y)

subplot(2,2,1), plot(x,y)   
subplot(2,2,3),plot(x,x-y)
hold on
subplot(2,2,1), plot(x-y,x)
% subplot will put multiple graphs 
% at a specified position in a grid

%%%%%%%%%%%%%
% miscellaneas 
%%%%%%%%%%%%%

num=1.23457809; 
format short %changes the number of digits shown 
num 
format long %changes the number of digits shown  
num  
%do help format   

clear all %clears all the variables 
close all %closes all the windows   

%%%%%%%%%%
% for loop 
%%%%%%%%%%  

i1=1; 
iN=5; 
step=2; 
for it=i1:step:iN 
    it^2   
end;      
    
%%%%%%%%%%%%%%  
% if & while statement  
%%%%%%%%%%%%%%   

i1=1; 
iN=5; 

if (i1<iN)  %try different conditions      
    i1=i1^2 
end;     

while(i1<iN)
    i1=i1^2
end;

%%%%%%%%%%%%%%  
% * and .*  
% note the difference between the two ways of doing multiplications 
% just * does the usual matrix matrix multiplication 
% .* does the multiplication of the corresponding elements 
%%%%%%%%%%%%%%   

J=[1,2;3,4]  
K=J*J 
L=J.*J 

f=[1,2,3,4]; 
g=f.*f %this is equivalent to a scalar product (if the vectors are real) 
h=f'*f %this is also equivalent to a scalar product but we take the adjoint of the vector first  



