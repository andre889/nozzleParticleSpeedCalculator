function [mp,xtemp] = NozzleMesh(x,y,div,ploting)
% Returns:  mp (y coordinates for each mesh point)
%           xtemp (x coordinates for each mesh point)
%
% Inputs:   x (x coordinates for each point on nozzle)
%           y (y coordinates for each point on nozzle)
%           div (divisions from zero to y at each point on nozzle)
%           plotting (choose to plot all points or not)
%
% Uses: For creating cells and meshing for simulating nozzle
%       performance in openfoam
% Objective: To create points at each vertex of a cell in a mesh to be used
%            in a future meshing function
% Autor: Devin McGee
% Project: MURI
% Note: written to be compatible with C++ (non vectorized)

for i=1:length(x)
    for j=1:div+1
        mp(j,i) = y(i)*(j-1)/div;
    end
end

for n=1:length(x)
    for i=1:div+1
        xtemp(i,n) = x(n);
    end
end

if strcmp(ploting,'Yes')==1||strcmp(ploting,'yes')==1
    hold on
    for k=1:length(x)
        for m=1:div+1
            plot(xtemp(m,k),mp(m,k),'r*','markersize',1)
        end
    end
    hold off
end
end