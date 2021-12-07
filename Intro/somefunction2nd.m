% this is a simple example of function of a first order ode    
function dydt=somefunction2nd(t,y) 
    dydt_1 = y(2); % here set the derivative of the first element of the vector y to be equal to the second element   
    dydt_2 = -y(1); % here you write the relevant d^2y/dt^2    here we chose a simple example which will give sinuisoidal evolution       
    dydt=[dydt_1; dydt_2];  
end  







