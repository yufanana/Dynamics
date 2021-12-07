%example of function with two inputs and two outputs 
% the function itself has a nested function

function [output1,output2]=function_name(input1, input2)  

output1=input1*input2; 
output2=innerfunction(output1); % calls another function 

    function out=innerfunction(in) % function defined within the function. It can also be defined outside 
    out=in.^2; 
    end 

end 

