function [x, A, A_star,ThroatLoc,y] = getNozzle6(leng,Ed,Post_Length,Throat_Radius, endlength, divpm, minpt)
% Returns:  x (position [m] along the nozzle)
%           A (Area [m^2] along the nozzle at each x)
%           A_Star (Area [m^2] of the throat)
%           ThroatLoc (Location of the throat [m])
%
% Inputs:   leng (Desired lenght of the throat [m])
%           Ed (exit radius [m])
%           Post_Length (Lenght of expansion section after throat
%           requirements [m])
%           Throat_Radius (Radius of the Troat [m])
%           endlength (straight section lenght after expansion [m])
%           divpm (divisions per meter in nozzle, helpful for mesh)
%
% For information on the nozzle geometry chosen (Rao) https://doi.org/10.1016/j.heliyon.2019.e01273
% For information on the equations required to draw the rao nozzle http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
% Author: Devin McGee 
% Project: MURI
% Use: First Draft Nozzle Built By Shop

thN = 22; % angle of nozzle after throat requirements
thE = 0; % angle of nozzle after expansion
Div1 = max(ceil(0.00026517*divpm),minpt);
Div2 = max(ceil((0.00053033-0.00026517)*divpm),minpt);
Div3 = max(ceil((0.00056611-0.00053033)*divpm),minpt);
th1 = linspace(-135,-90, Div2); % Set the bounds for first curve to be drawn
th2 = linspace(-90,(thN-90), Div3); % Set the bounds for the second curve to be drawn
Rt = Throat_Radius; % set the throat radius to used variable
A_star =  pi*(Rt)^2; % calculate the throat area
xin = 1.5*Rt*cosd(th1); % calculate the x position of the first nozzle curve (eq 4 in nozzle design)
xin = xin+abs(xin(1)); % move the start point to zero
yin = 1.5*Rt*sind(th1)+Rt*1.5+Rt;% calculate the y position of the first nozzle curve (eq 4 in nozzle design)

th3 = linspace(45,90, Div1); % set the bounds for an added inlet section
xinp = 1.5*Rt*cosd(th3)+0.25*xin(1); % calculate the extra inlet x positions
yinp = 1.5*Rt*sind(th3)+yin(1); % calculate the extra inlet x positions
xinp = 1.5*Rt*cosd(th3)+0.25*xin(1)+(xin(1)-xinp(1)); % move the x to match xin start point
yinp = 1.5*Rt*sind(th3)+yin(1)-(yinp(1)-yin(1)); % move the y to match yin start point
xinp = flip(xinp); % make the points asscend instead of descend
yinp = flip(yinp); % make the points asscend instead of descend

xin = xin+abs(xinp(1)); % move the xin positions to the new zero
xinp = xinp+abs(xinp(1)); % move the xinp to the new xin start
ThroatLoc = xin(end); % find the throat location
xse = xin(end)+leng; % find the position of the end of the nozzle straight section
if leng==0
    xs = xin(end);
    ys = yin(end);
else
    xs = linspace(xin(end),xse,max(ceil(leng*divpm),minpt)); % make the nozzle straight section x positions
    ys = ones(size(xs))*yin(end); % make the nozzle straight section y positions
end
y2 = 0.382*Rt*sind(th2)+0.382*Rt+Rt; % Make the initial expansion section y poisitions (equaiton 5)
x2 = 0.382*Rt*cosd(th2); % Make the initial expansion section x poisitions (equaiton 5)
x2 = x2+abs(xs(end)); % Set the start point to be at the end of the straight section
plot(x2,y2); % plot the nozzle
hold on; % plot the nozzle
plot(xin,yin); % plot the nozzle
plot(xinp,yinp); % plot the nozzle
plot(xs,ys); % plot the nozzle
plot(ThroatLoc,yin(end),'r*') % plot the throat location

Nx = x2(end); % Set up for the expansion section calculation (equations 6-10)
Ny = y2(end); % Set up for the expansion section calculation (equations 6-10)
Ex = leng+Post_Length; % Set up for the expansion section calculation (equations 6-10)
Ey = Ed; % Set up for the expansion section calculation (equations 6-10)
m1 = tand(thN); % Set up for the expansion section calculation (equations 6-10)
m2 = tand(thE); % Set up for the expansion section calculation (equations 6-10)
c1 = Ny-Nx*m1; % Set up for the expansion section calculation (equations 6-10)
c2 = Ey-m2*Ex; % Set up for the expansion section calculation (equations 6-10)
Qx = (c2-c1)/(m1-m2); % Set up for the expansion section calculation (equations 6-10)
Qy = (m1*c2-m2*c1)/(m1-m2); % Set up for the expansion section calculation (equations 6-10)

t = linspace(0,1,max(ceil(Post_Length*divpm),minpt)); % Set up for the expansion section calculation (equations 6-10)

xbell = ((1-t).^2.*Nx + 2*(1-t).*t.*Qx + t.^2.*Ex); % Calulate the expansion section x position (equation 6)
ybell = ((1-t).^2.*Ny + 2*(1-t).*t.*Qy + t.^2.*Ey); % Calulate the expansion section y position (equation 6)

plot(xbell,ybell) % plot the nozzle

% Ar = (pi*Ey^2)/(pi*Rt^2);
if endlength==0
    xend = [];
    yend = [];
else
    xend = linspace(xbell(end),(xbell(end)+endlength),max(ceil(endlength*divpm),minpt)); % add a straight section to the end
    yend = ones(size(xend))*ybell(end); % add a straight section to the end
end
plot(xend,yend) % plot the straigh section
x = [xinp, xin(2:end), x2(2:end), xbell(2:end), xend(2:end)]; % combine all the x locations for exporting
y = [yinp, yin(2:end), y2(2:end), ybell(2:end), yend(2:end)]; % combine all the y locations for exporting
A = pi*y.^2; % Calculate the area at all y locations 
hold off

    
end