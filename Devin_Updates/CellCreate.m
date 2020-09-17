function [cell] = CellCreate(xpos,ypos,thick)
% Returns:  cell (structure containing mesh face properties)
%               cell.vert (verticies of each face in x,y,z)
%               cell. 
%
% Inputs:   xpos (x coordinates for each point on nozzle)
%           ypos (y coordinates for each point on nozzle)
%
% Uses: For creating cells and meshing for simulating nozzle
%       performance in openfoam
% Objective: To define all properties required by openfoam for defining the
%            faces of mesh 
% Autor: Devin McGee
% Project: MURI
% Note: written to be compatible with C++ (non vectorized)

clear cell
% for i=1:length(xpos(1,:))
%     for j=1:length(xpos(:,1))
%         cell.point{pos,1} = [xpos(1,i), ypos(1,j),0]; % start at begining and count forward
%         cell.point{pos,1} = [xpos(1,i), ypos(1,j),0]; % start at end and count back
%         pos = pos+1;
%     end
% end

pos = 1;
for i=1:length(xpos(:,1))
    cell.point{pos,1} = [xpos(i,1), ypos(i,1), 0];
    pos = pos + 1;
    cell.point{pos,1} = [xpos(i,1), ypos(i,1), thick];
    pos = pos + 1;
end


for j=2:length(xpos(1,:))
    cell.point{pos,1} = [xpos(1,j), ypos(end,j), 0];
    pos = pos + 1;
    cell.point{pos,1} = [xpos(1,j), ypos(end,j), thick];
    pos = pos + 1;
end

ypos(:,end) = flip(ypos(:,end));
for i=2:(length(xpos(:,1)))
    cell.point{pos,1} = [xpos(i,end), ypos(i,end), 0];
    pos = pos + 1;
    cell.point{pos,1} = [xpos(i,end), ypos(i,end), thick];
    pos = pos + 1;   
end

xpos(1,:) = flip(xpos(1,:));
ypos(:,end) = flip(ypos(:,end));
for j=2:(length(xpos(1,:))-1)
    cell.point{pos,1} = [xpos(1,j), ypos(1,j), 0];
    pos = pos + 1;
    cell.point{pos,1} = [xpos(1,j), ypos(1,j), thick];
    pos = pos + 1;
end
xpos(1,:) = flip(xpos(1,:));

% pos = 1;
% for i=1:(length(xpos(1,:))-1)
%     for j=1:(length(xpos(:,1))-1)
%         cell.vert{pos,1} = [xpos(1,i+1),ypos(j,i),0];
%         cell.vert{pos,2} = [xpos(1,1+1),ypos(j+1,i),0];
%         cell.vert{pos,3} = [xpos(1,i),ypos(j+1,i),0];
%         cell.vert{pos,4} = [xpos(1,i),ypos(j,i),0];
%         cell.vert{pos,1} = [xpos(1,i+1),ypos(j,i),thick];
%         cell.vert{pos,2} = [xpos(1,1+1),ypos(j+1,i),thick];
%         cell.vert{pos,3} = [xpos(1,i),ypos(j+1,i),thick];
%         cell.vert{pos,4} = [xpos(1,i),ypos(j,i),thick];
%         
%         if j==1||j==(length(xpos(:,1))-1)
%             cell.bound{pos,1} = 'Boundry Cell';
%         else
%             cell.bound{pos,1} = 'Face Cell';
%         end
%         
%         pos = pos+1;
%         
%     end
% end
% 
% post = pos-1;
% for k=1:4
%     for m=1:post
%         cell.face(m,k) = m + post*(k-1);
%     end
% end

end


