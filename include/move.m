%{
Updates the x- and y-coordinates of a residue subjected to a trial move.
Returns the updated x- and y-coordinates of all residues.

IN:  vector representing x-coordinate of all residues;
     vector representing y-coordinate of all residues;
     index of residue subjected to a trial move;
     magnitude of trial movement of residue in x- and y-axes.
OUT: updated vector representing x-coordinate of all residues;
     updated vector representing y-coordinate of all residues.
%}
function [xvals, yvals] = move(xvals, yvals, trial_pos, x_moved, y_moved)
    xvals(trial_pos) = xvals(trial_pos) + x_moved;
    yvals(trial_pos) = yvals(trial_pos) + y_moved;
end