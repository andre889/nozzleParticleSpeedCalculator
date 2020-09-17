%Return x and A of nozzle design 1
%currently a test case
%Author: Austin Andrews, 7/27/20
%Project: MURI
%http://www.aerorocket.com/Nozzle/Validate/Validate.html

function [x, A, A_star,ThroatLoc] = getNozzle5()
    resolution = 400; 
    x = linspace(0,20e-2,resolution); % in [m]
    A_star =  pi*(0.2e-3)^2/4; %throat size in [m]
    A = x*0.0;
    ThroatLoc = 25e-3;
    Ainlet = 4.4e-3^2*pi/4;
    Aoutlet = 2.0e-3^2*pi/4;
    A1 = (Ainlet-A_star)/((ThroatLoc^2));
    A2 = (Aoutlet-A_star)/((max(x)-ThroatLoc)^2);
    A(x > ThroatLoc) = A2*(x(x > ThroatLoc)-ThroatLoc).^2 + A_star; 
    A(or(x <ThroatLoc,x == ThroatLoc)) = A1*(x(or(x <ThroatLoc,x == ThroatLoc))-ThroatLoc).^2 + A_star;
    %debug
    figure
    plot(x,A)
    xlabel('Distance [m]')
    ylabel('Area [m^2]')
end 