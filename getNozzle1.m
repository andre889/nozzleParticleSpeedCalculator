%Return x and A of nozzle design 1
%currently a test case
%Author: Austin Andrews, 7/27/20
%Project: MURI
%http://www.aerorocket.com/Nozzle/Validate/Validate.html

function [x, A, A_star,ThroatLoc] = getNozzle1()
    resolution = 100; 
    x = linspace(0,1.0,resolution); % in [m]
    A_star =  .001387; %throat size in [m]
    A = x*0;
    ThroatLoc = 0.25;
    A(x > ThroatLoc) = 0.1089*(x(x > ThroatLoc)-ThroatLoc).^2 + A_star; 
    A(or(x <ThroatLoc,x == ThroatLoc)) = 0.1089*(x(or(x <ThroatLoc,x == ThroatLoc))-ThroatLoc).^2 + A_star;
    %debug
    figure
    plot(x,A)
    xlabel('Distance [m]')
    ylabel('Area [m^2]')
end 