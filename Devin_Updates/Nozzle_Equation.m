%% Nozzle shape equation

thN = 22;
thE = 0;
th1 = linspace(-135,-90);
th2 = linspace(-90,(thN-90));
Rt = 0.25;
xin = 1.5*Rt*cosd(th1);
yin = 1.5*Rt*sind(th1)+Rt*1.5+Rt;
plot(xin,yin)

th3 = linspace(45,90);
xinp = 1.5*Rt*cosd(th3)+0.25*xin(1);
yinp = 1.5*Rt*sind(th3)+yin(1);
xinp = 1.5*Rt*cosd(th3)+0.25*xin(1)+(xin(1)-xinp(1));
yinp = 1.5*Rt*sind(th3)+yin(1)-(yinp(1)-yin(1));
y2 = 0.382*Rt*sind(th2)+0.382*Rt+Rt;
x2 = 0.382*Rt*cosd(th2);
plot(x2,y2); hold on; plot(xin,yin); plot(xinp,yinp)

Nx = x2(end);
Ny = y2(end);

Ex = 300;
Ey = 2;


m1 = tand(thN);
m2 = tand(thE);
c1 = Ny-Nx*m1;
c2 = Ey-m2*Ex;
Qx = (c2-c1)/(m1-m2);
Qy = (m1*c2-m2*c1)/(m1-m2);

t = linspace(0,1,1000);

xbell = ((1-t).^2.*Nx + 2*(1-t).*t.*Qx + t.^2.*Ex);
ybell = ((1-t).^2.*Ny + 2*(1-t).*t.*Qy + t.^2.*Ey);

% plot(xbell,ybell)

Ar = (pi*Ey^2)/(pi*Rt^2);

% axis([-1 300 0 2])

hold off



