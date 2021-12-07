% this function contains the ode to be solved     	
function dXdt = W4_CP7_slope_dxdt(T,X,g,k,A,m)
    W4_CP7_sinacosa=W4_CP7_angle(X(1)); 
        dXdt = zeros(2,1); 
        dXdt = [X(2); (-k/m*X(2)+g*W4_CP7_sinacosa)/(1+A)];   
    end