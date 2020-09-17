%% Test case for meshing in openFoam
% Written to test meshing program with exisiting openFoam software. This
% test case will be used to mesh nozzles in the future. 

%% General Setup
clear
Filep = 'Meshing_Output_Points.csv';
Filef = 'Meshing_Output_Faces.csv';
Fileo = 'Meshing_Output_Owner.csv';
Filen = 'Meshing_Output_Neighbor.csv';
VertStep = 10;
HorStep = 10;
Depth = 1e-4;
nFace = 4; % number of points on a Face (really only works with 4)

a = [1 2 3; 4 5 6; 7 8 9]; %... because

%% Create the geometry
% Create the 2D face on the origin
Leng = 0.05; % length in meters
Height = 0.01; % height in meters

xpos = linspace(0, Leng, HorStep);
ypos = zeros(VertStep, length(xpos));
for i=1:length(xpos)
    ypos(:,i) = linspace(0, Height, VertStep);
end

xposa = zeros(VertStep, length(xpos));
for i=1:length(ypos)
    xposa(:,i) = xpos(i)*ones(VertStep,1);
end

%% Import the nozzle geometery
% xtemp = load('nozzle_xpos');
% ytemp = load('nozzle_ypos');
% 
% xposa = xtemp.xtemp;
% ypos = ytemp.mp;
%% Label the points
pos = 1;
for i=1:length(xposa(:,1))
    cell.point(pos,:) = [xposa(i,1), ypos(i,1), 0];
    pos = pos + 1;
    cell.point(pos,:) = [xposa(i,1), ypos(i,1), Depth];
    pos = pos + 1;
end


for j=2:length(xposa(1,:))
    cell.point(pos,:) = [xposa(1,j), ypos(end,j), 0];
    pos = pos + 1;
    cell.point(pos,:) = [xposa(1,j), ypos(end,j), Depth];
    pos = pos + 1;
end

ypos(:,end) = flip(ypos(:,end));
for i=2:(length(xposa(:,1)))
    cell.point(pos,:) = [xposa(i,end), ypos(i,end), 0];
    pos = pos + 1;
    cell.point(pos,:) = [xposa(i,end), ypos(i,end), Depth];
    pos = pos + 1;   
end

xposa(1,:) = flip(xposa(1,:));
ypos(:,end) = flip(ypos(:,end));
for j=2:(length(xposa(1,:))-1)
    cell.point(pos,:) = [xposa(1,j), ypos(1,j), 0];
    pos = pos + 1;
    cell.point(pos,:) = [xposa(1,j), ypos(1,j), Depth];
    pos = pos + 1;
end
xposa(1,:) = flip(xposa(1,:));

for i=2:(length(xposa(1,:))-1)
    if rem(i,2)==0
        for j=2:(length(ypos(:,1))-1)
            cell.point(pos,:) = [xposa(j,i), ypos(j,i), 0];
            pos = pos + 1;
        end
    else
        ypos(:,i) = flip(ypos(:,i));
        for j=2:(length(ypos(:,1))-1)
            cell.point(pos,:) = [xposa(j,i), ypos(j,i), 0];
            pos = pos + 1;
        end
        ypos(:,i) = flip(ypos(:,i));
    end
end

for i=2:(length(xposa(1,:))-1)
    if rem(i,2)==1
        for j=2:(length(ypos(:,1))-1)
            cell.point(pos,:) = [xposa(j,i), ypos(j,i), Depth];
            pos = pos + 1;
        end
    else
        ypos(:,i) = flip(ypos(:,i));
        for j=2:(length(ypos(:,1))-1)
            cell.point(pos,:) = [xposa(j,i), ypos(j,i), Depth];
            pos = pos + 1;
        end
        ypos(:,i) = flip(ypos(:,i));
    end
end

%% Create the Faces and Find the face centers for use in ownership
% Lets go Vertical First

faces = 1;
for j=1:length(xposa(1,:))
    for i=1:(length(xposa(:,1))-1)
        if i==0
            % do litterally nothing.. not great code but yeah
        else
            p1temp = [xposa(1,j), ypos(i,1), 0];
            p2temp = [xposa(1,j), ypos(i,1), Depth];
            p3temp = [xposa(1,j), ypos(i+1,1), 0];
            p4temp = [xposa(1,j), ypos(i+1,1), Depth];
            for m=1:length(cell.point(:,1))
                [~,r1] = ismember(p1temp,cell.point,'rows');
                [~,r2] = ismember(p2temp,cell.point,'rows');
                [~,r3] = ismember(p3temp,cell.point,'rows');
                [~,r4] = ismember(p4temp,cell.point,'rows');
            end
            cell.face(faces,:) = [r1, r2, r3, r4];
            cell.fcent(faces,:) = [xposa(1,j), (ypos(i,1)+ypos(i+1,1))/2, (0+Depth)/2];
            faces = faces + 1;
        end
    end
end

flaghor = faces;
% Then Lets go Horizontal
for j=1:length(xposa(:,1))
    for i=1:(length(xposa(1,:))-1)
        if i==0
            % do litterally nothing.. not great code but yeah
        else
            p1temp = [xposa(1,i), ypos(j,1), 0];
            p2temp = [xposa(1,i), ypos(j,1), Depth];
            p3temp = [xposa(1,i+1), ypos(j,1), 0];
            p4temp = [xposa(1,i+1), ypos(j,1), Depth];
            for m=1:length(cell.point(:,1))
                [~,r1] = ismember(p1temp,cell.point,'rows');
                [~,r2] = ismember(p2temp,cell.point,'rows');
                [~,r3] = ismember(p3temp,cell.point,'rows');
                [~,r4] = ismember(p4temp,cell.point,'rows');
            end
            cell.face(faces,:) = [r1, r2, r3, r4];
            cell.fcent(faces,:) = [(xposa(1,i)+xposa(1,i+1))/2, ypos(j,1), (0+Depth)/2];
            faces = faces + 1;
        end
    end
end

flagback = faces;
% Alright now I Guess the Front and Back Planes
for j=1:(length(xposa(:,1))-1)
    for i=1:(length(xposa(1,:))-1)
        if i==0
            % do litterally nothing.. not great code but yeah
        else
            p1temp = [xposa(1,i), ypos(j,1), 0];
            p2temp = [xposa(1,i+1), ypos(j,1), 0];
            p3temp = [xposa(1,i+1), ypos(j+1,1), 0];
            p4temp = [xposa(1,i), ypos(j+1,1), 0];
            for m=1:length(cell.point(:,1))
                [~,r1] = ismember(p1temp,cell.point,'rows');
                [~,r2] = ismember(p2temp,cell.point,'rows');
                [~,r3] = ismember(p3temp,cell.point,'rows');
                [~,r4] = ismember(p4temp,cell.point,'rows');
            end
            cell.face(faces,:) = [r1, r2, r3, r4];
            cell.fcent(faces,:) = [(xposa(1,i)+xposa(1,i+1))/2, (ypos(j,1)+ypos(j+1,1))/2, 0];
            faces = faces + 1;
        end
    end
end

flagfront = faces;
for j=1:(length(xposa(:,1))-1)
    for i=1:(length(xposa(1,:))-1)
        if i==0
            % do litterally nothing.. not great code but yeah
        else
            p1temp = [xposa(1,i), ypos(j,1), Depth];
            p2temp = [xposa(1,i+1), ypos(j,1), Depth];
            p3temp = [xposa(1,i+1), ypos(j+1,1), Depth];
            p4temp = [xposa(1,i), ypos(j+1,1), Depth];
            for m=1:length(cell.point(:,1))
                [~,r1] = ismember(p1temp,cell.point,'rows');
                [~,r2] = ismember(p2temp,cell.point,'rows');
                [~,r3] = ismember(p3temp,cell.point,'rows');
                [~,r4] = ismember(p4temp,cell.point,'rows');
            end
            cell.face(faces,:) = [r1, r2, r3, r4];
            cell.fcent(faces,:) = [(xposa(1,i)+xposa(1,i+1))/2, (ypos(j,1)+ypos(j+1,1))/2, Depth];
            faces = faces + 1;
        end
    end
end

% make a list of face names
for i=1:faces-1
    face(i,1) = i;
end

%% Give ownership
% thoughts: give faces a center point, give cells a center point
% use cell center point to find 6 closest faces and label them a owned
% use cell center to find the neighbor cells

cent = 1;
for i=1:(length(xposa(:,1))-1)
    for j=1:(length(xposa(1,:))-1)
        cell.cent(cent,:) = [(xposa(j,i)+xposa(j,i+1))/2, (ypos(j,i)+ypos(j+1,i))/2, (0+Depth)/2];
        cent = cent + 1;
    end
end

for i=1:length(cell.cent(:,1))
    distance{i,1} = [abs(cell.cent(i,1)-cell.fcent(:,1)), abs(cell.cent(i,2)-cell.fcent(:,2)), abs(cell.cent(i,3)-cell.fcent(:,3))];
end

% create the first cell (cell at origin)
ownv = owner(1, flaghor, distance{1,1}, 2);
outh = owner(flaghor, flagback, distance{1,1}, 2);
outfr = owner(flagback, flagfront, distance{1,1}, 2);
outb = owner(flagfront, length(distance{1,1})+1, distance{1,1}, 2);

cell.own{1,1} = [ownv outh outfr outb];

% create the first row of cells (vertical from origin), lose one horizontal
% face of ownership
for i=2:(length(xposa(:,1))-1)
    ownv = owner(1, flaghor, distance{i,1}, 2);
    outh = owner(flaghor, flagback, distance{i,1}, 1);
    outfr = owner(flagback, flagfront, distance{i,1}, 2);
    outb = owner(flagfront, length(distance{i,1})+1, distance{i,1}, 2);
    cell.own{i,1} = [ownv outh outfr outb];
end
endpoint = i-1; % subtract one so the next loop can start at 2 but still but a point after the last i

% Create next row of cells (horizontal from origin), lose one vertical face
% of ownership
for i=2:(length(xposa(1,:))-1)
    j = (i-1)*(length(xposa(:,1))-1)+1;
    ownv = owner(1, flaghor, distance{j,1}, 1);
    outh = owner(flaghor, flagback, distance{j,1}, 2);
    outfr = owner(flagback, flagfront, distance{j,1}, 2);
    outb = owner(flagfront, length(distance{j,1})+1, distance{j,1}, 2);
    cell.own{i+endpoint,1} = [ownv outh outfr outb];
end
endpoint = endpoint+i; 

% Now create the rest of the cells, all cells will lose one vertical and
% one horizontal face of ownership
count = 1;
for i=2:(length(xposa(1,:))-1)
    for j=2:(length(xposa(:,1))-1)
        k = length(xposa(:,1))*(i-1)+1+(j-2)-(i-2); % this could be a problem area for big mesh.. need to check
        temp(j,i) = k;
        ownv = owner(1, flaghor, distance{k,1}, 1);
        outh = owner(flaghor, flagback, distance{k,1}, 1);
        outfr = owner(flagback, flagfront, distance{k,1}, 2);
        outb = owner(flagfront, length(distance{k,1})+1, distance{k,1}, 2);
        cell.own{count+endpoint,1} = [ownv outh outfr outb];
        count = count + 1;
    end
end

%% Tell the cells what the neighbors are

% This makes the assumption all cells are on the same plane
% A tolerance will be used much like the owner funciton to relieve any
% matlab rounding issues
tol = 1e-10;
% i = 10;

for i=1:length(cell.cent(:,1))
    % this finds the x and y distance between the cell of interest and all
    % the cells in the domain
    dif1 = abs(cell.cent(:,1)-cell.cent(i,1));
    dif2 = abs(cell.cent(:,2)-cell.cent(i,2));
    % this finds the minimum distance for x and y in the difference. The
    % cells I want are the ones in the x and y row 
    min1 = min(dif1);
    min2 = min(dif2);
    % This tells me all the cells in the x and y row and column I am going to be working
    % with 
    col1 = find(abs(min1-dif1)<tol);
    col2 = find(abs(min2-dif2)<tol);
    
    % the goal here is to use the x and y row and column of cells and use
    % my current cell to say which of the cells in the x and y row and
    % column are the colsest and which are not. The challenge was cells
    % having different number of neighbors and that led to having lots of
    % if else statements. 
    if find(col1==i)==1&&find(col2==i)>1&&i~=col2(end)
        % this makes the bottom face of neighbors and ignores the x=end
        % neighbors
        yneigh = col1(find(col1==i)+1);
        xneigh = [col2(find(col2==i)-1) col2(find(col2==i)+1)];
    elseif i==col1(end)&&find(col2==i)==1
        % this makes the left side top corner cell
        xneigh = col2(find(col2==i)+1);
        yneigh = col1(find(col1==i)-1);
    elseif i==col1(end)&&i==col2(end)
        % this makes the right side top corner cell
        xneigh = col2(find(col2==i)-1);
        yneigh = col1(find(col1==i)-1);
    elseif i==col1(end)&&find(col2==i)>1
        % this makes the top row of cells
        yneigh = col1(find(col1==i)-1);
        xneigh = [col2(find(col2==i)-1) col2(find(col2==i)+1)];
    elseif find(col2==i)==1&&find(col1==i)>1
        % this makes the left face of cells 
        xneigh = col2(find(col2==i)+1);
        yneigh = [col1(find(col1==i)-1) col1(find(col1==i)+1)];
    elseif find(col2==i)==1&&find(col1==i)==1
        % this makes the very first cell (origin cell)
        xneigh = col2(find(col2==i)+1);
        yneigh = col1(find(col1==i)+1);
    elseif i==col2(end)&&find(col1==i)==1
        % this makes the bottom right corner cell)
        xneigh = col2(find(col2==i)-1);
        yneigh = col1(find(col1==i)+1);
    elseif i==col2(end)&&find(col1==i)>1
        % this makes the right wall of cells 
        xneigh = col2(find(col2==i)-1);
        yneigh = [col1(find(col1==i)-1) col1(find(col1==i)+1)];
    else
        % this makes all the middle cells 
        yneigh = [col1(find(col1==i)-1) col1(find(col1==i)+1)];
        xneigh = [col2(find(col2==i)-1) col2(find(col2==i)+1)];
    end
    % this writes to the variable I want at the very end
    cell.neigh{i,1} = [xneigh yneigh];
end



%% File output

% append face numbers to face file
cell.face = [face, cell.face];

% Write to file
writematrix(cell.point,Filep)
writematrix(cell.face, Filef)
writecell(cell.own, Fileo)
writecell(cell.neigh, Filen)












