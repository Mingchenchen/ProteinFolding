%{
Defines distance for moving any residue in terms of x- and y-axes for a
randomly chosen direction between 1 and 4 (inclusive).
Let 1 be N, 2 be S, 3 be E, 4 be W.

OUT: distance to be moved by some residue along x-axis;
     distance to be moved by some residue along x-axis.
%}
function [x_moved, y_moved] = getDistancesToBeMoved()
    direction = randi(4);

    NORTH = 1;
    SOUTH = 2;
    EAST  = 3;
    WEST  = 4;
    
    if direction == NORTH
        x_moved = 0; y_moved = 1;
    elseif direction == SOUTH
        x_moved = 0; y_moved = -1;
    elseif direction == EAST
        x_moved = 1; y_moved = 0;
    elseif direction == WEST
        x_moved = -1; y_moved = 0;
    end
end