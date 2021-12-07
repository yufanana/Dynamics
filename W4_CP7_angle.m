% this function returns sin(alpha)cos(alpha)    
function W4_CP7_sinacosa = W4_CP7_angle(x)
    eps=0.0001; 
    tana=-(W4_CP7_traj(x+eps)-W4_CP7_traj(x))/eps; 
    W4_CP7_sinacosa= tana./(1+tana.^2); 
 