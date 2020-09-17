%This program will be used to calculate the gas properties in a supersonic
%nozzle using isentropic relationships as found here 
%https://www.grc.nasa.gov/WWW/k-12/airplane/nozzled.html
%Using these properties, particle impact velocites are determined using the
%drag expression from Li et al. https://doi.org/10.1016/j.ijheatmasstransfer.2018.10.028
%Author: Austin Andrews, 7/27/20
%Project: MURI

close all
clear all
%particle Properties
rhoP = 1000; %density of particle [kg/m^3]
dp = exp(linspace(log(1e-7),log(10e-6),100)); 

%Gas properties
% gamma = 1.4; %Specific Heat ratio
gamma = 1.667; %Specific Heat ratio helium
% M = 28.9/1000; %Gas Molecular Weight [kg/Mol]
M = 4.0/1000; %helium Gas Molecular Weight [kg/Mol]
R = 8.314/M; % gas Specific R [J/(kg*K)]
Tt = 500; %total temperature [K] *at inlet
Pt = 101325*.9; %total pressure [pa] *at inlet
rhoT = Pt/(R*Tt); %total density [kg/m^3] *at inlet
%using A/A_star = ((gamma+1)/2)^(-(gamma+1)/2/(gamma-1)) * (1 + 
%M^2(gamma-1)/2)^((gamma+1)/(2*(gamma-1))/M 

%given a function That 1. determines the distance x[m] and 2. determins
%A[m^2] as a function of x find the Mach number M
leng = 0.0; % Throat Lenght in m
Expansion = 2.0;
Throat_Radius = 1*1e-3; % Throat Radium in m
Ed = Throat_Radius*Expansion; % exit radius in m
Post_Length = 40e-3; % lenght of nozzle after throat in m
endlength = 80e-3; % length of the nozzle after expansion
divpm = 250; % divisisions per meter
minpt = 15;
taper = 10;
inlet = 100e-3;
[x, A, A_star,ThroatLoc,y] = getNozzle7(leng,Ed,Post_Length, Throat_Radius, endlength,divpm,minpt, taper, inlet);
% div = 29; % divisions along y axis for each x for mesh (note, extra point added at 0 so div+1 is result)
% ploting = 'no';
% [mp,xtemp] = NozzleMesh(x,y,div,ploting);

% CellThickNess = 1e-5;
% [cell] = CellCreate(xtemp,mp,CellThickNess);
%newton method to solve for Mach number 
[Msub,Msup] = sub_super(A./A_star,gamma);

%Stitch together the two solutions 
Ma(x >= ThroatLoc) = Msup(x >= ThroatLoc);
Ma(x < ThroatLoc) = Msub(x < ThroatLoc);

%Solve for other nozzle variables 
T = Tt*(1 + 0.5*(gamma-1).*Ma.*Ma).^(-1.0); % Temperature at x position
P = Pt*(1 + 0.5*(gamma-1).*Ma.*Ma).^(-gamma/(gamma-1)); % Pressure at x position
%mfp = getMFP(T,P); %Assumed air 
a = sqrt(gamma*(8.314/M)*T);  %speed of sound 
U = Ma.*a; % Gas Velocity 
rho = P./((8.314/M).*T);
mdot = A_star*sqrt(gamma/R)*((gamma+1)/2)^(-1*(gamma+1)/2/(gamma-1))*Pt/sqrt(Tt) %mass flow in kg/s;
% mdot = A_star*sqrt(gamma*(2/(gamma+1))^((gamma+1)/(gamma-1))*Pt*rhoT); %mass flow in kg/s;
% mdot = Pt*A_star*sqrt(gamma/(R*Tt))*(1+(gamma-2)/2)^((gamma+1)/(-2*(gamma-1)));
Vcylinder = 6.14476; % volume of Helium in Cylinder
Water_Vol = 0.0499; % Volume of water in cylinder
Pg = Pt*Vcylinder/Water_Vol; % Pressure of Helium in Cylinder
CylMass = Pg/(273*R)*Water_Vol; % Mass of Helium in Cylinder
Vdot = mdot/rho(1); % m^3/s
VdotLPM = 60000*Vdot
rhoSTP = 101325./(R.*273); % Desity of Gas at STP
SLPM = 60000*mdot/rhoSTP % Flow rate of Gas at STP
% Cylinders_per_Hour = (Vcylinder/Vdot)/60/60
Cylinders_per_Hour = mdot/CylMass/60/60 % Cylinders of Gas used per hour
%mu = 1.82e-5 * ((273 + 110.4)./(T+110.4)).* (T./273).^(3/2); %viscosity sutherland need ref also for other gases?
mu = 1.87e-5*(T/273).^0.668; %power law helium see https://en.wikipedia.org/wiki/Temperature_dependence_of_viscosity
mfp = (mu./P).*sqrt(pi*8.314.*T./2.0/M); %see https://en.wikipedia.org/wiki/Mean_free_path

%% particles 

particle_velocity = dp*0;
for i = 1:length(dp)
Cc = getCc(dp(i),mfp);
y0 = [0 U(1)]; %initial value same as gass velocity
tspan = [0 0.01];
%options = odeset('RelTol',1e-8,'AbsTol',1e-10);
%Opt    = odeset('Events', @myEvent);
opts = odeset('RelTol',1e-4,'AbsTol',1e-4,'Events', @myEvent,'Stats','on','InitialStep',1e-8);
[t,y] = ode45(@(t,y) odefcn(t,y,x,U,mu,mfp,dp(i)),tspan,y0,opts);

particle_velocity(i) = y(end,2); 
i
end

%% Plotting

figure
plot(dp*1e9,particle_velocity/U(end))
xlabel('Paritcle Diameter [nm]')
ylabel('Particle Velocity / Exit Velocity')

figure
plot(x,Ma)
xlabel('Distance [m]')
ylabel('Ma')

figure
plot(x,T)
xlabel('Distance [m]')
ylabel('Temperature [K]')

figure 
plot(x,U)
xlabel('Distance [m]')
ylabel('U [m/s]')


function [value, isterminal, direction] = myEvent(t, y)
value = (y(1) < 1); %nozzle length
isterminal = 1;   % Stop the integration
direction  = 0;
end

function dydt = odefcn(t,y,x,u,mu,mfp,dp);
 %y(1) is position
 %y(2) is velocity
 %u is gas velocity
 
 %Fd is drag force
 %Fd = 18*mu*Cd/(rhoP*dp*dp);
 mu_interp = interp1(x,mu,y(1),'linear','extrap');
 mfp_interp = interp1(x,mfp,y(1),'linear','extrap');
 
 rhoP = 1000;
 
 %Fd = U*18.0*visc*Cd*Re/Cc/rho/dp/dp/24.0;
 %Cd = 24/Re ; 
 Cc = getCc(dp,mfp_interp);
 Fd = 18.0*mu_interp/Cc/rhoP/dp/dp;
 
 dydt = zeros(2,1);
 Uint = interp1(x,u,y(1),'linear','extrap');
 dydt(1) = y(2);
 accel = -(y(2)-Uint)*Fd;
 dydt(2) = accel;
end
