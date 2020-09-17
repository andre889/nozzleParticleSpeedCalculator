function [out] = owner(start, flag, distance, f2take)
% Funciton to find which face centers are closest to the center of a cell
% of interest. 
% Inputs:
%        start: start of the range of interst (face number) thought is to
%           only pay attention to faces oriented the correct direction. i.e if
%           I want to find the vertical faces owned by a cell, the horizontal
%           faces could mess up my results by being closer to the center. 
%        flag: end of range of interst. see start for reasoning for this. 
%        distance: array of distance from cell of interest to all faces in
%           domain with x, y, z distance. 
%        f2take: number of faces to take, important for knowing how many
%           vertical and horizontal faces with closest centers actually belong
%           to the cell and what does not.
%
% Outputs:
%        out: face numbers correlating to the faces blonging to a cell in
%           the specified direction i.e. vertical horizontal ect. (either
%           one or two faces)
%
% Author: Devin McGee
% Project: MURI Noozzle Meshing
% Purpose: Face ownership by cell
% Date of last Revision: 9/13/2020

% set a difference tolerance. This needs to be done to take care of any
% matlab rounding errors that could lead to two numbers that are equal
% being counted as not equal. 
tol = 1e-10;

% This finds the minimum in a set range of cells of interest for the x, y,
% z coordinates. 
min1 = min(distance(start:(flag-1),1));
min2 = min(distance(start:(flag-1),2));
min3 = min(distance(start:(flag-1),3));
col1 = find(abs(min1-distance(start:(flag-1), 1))<tol);
col2 = find(abs(min2-distance(start:(flag-1), 2))<tol);
col3 = find(abs(min3-distance(start:(flag-1), 3))<tol);

if length(col1)>length(col2)
    eq1 = find(ismember(col1,col2)==1)+col1(1)-1;
else 
    eq1 = find(ismember(col2,col1)==1)+col2(1)-1;
end

if length(eq1)>length(col3)
    eqf = find(ismember(eq1,col3)==1)'+eq1(1)-1;
else
    eqf = find(ismember(col3,eq1)==1)'+col3(1)-1;
end

if f2take==2
    out = eqf+start-1;
elseif f2take==1
    out = eqf+start-1;
    out(:,1) = [];
else
    error('Unknown command for faces owned by cell')
end

end