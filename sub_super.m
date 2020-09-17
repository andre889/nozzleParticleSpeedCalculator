%Thanks to David Hill https://www.mathworks.com/matlabcentral/answers/482921-mach-number-area-relation-several-inputs
%Modified by Austin Andrews 7/227/20

function [Msub,Msup] = sub_super(ARatio,g)%ARatio needs to be >1
%g   = 1.4;
gm1 = g-1;
gp1 = g+1;
Msub=zeros(size(ARatio));
Msup=zeros(size(ARatio));
for i=1:length(ARatio)
    problem.objective = @(M) (1/M^2)*(((2+gm1*M^2)/gp1)^(gp1/gm1))-ARatio(i)^2;    
    problem.solver    = 'fzero';                                              
    problem.options   = optimset(@fzero);                                      
    problem.x0 = [1e-6 1];                                                     
    Msub(i)       = fzero(problem);                                               
    problem.x0 = [1+1e-6 50];                                                   
    Msup(i)       = fzero(problem);  
end

end

