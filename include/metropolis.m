%{
Metropolis constraint: if helper_exponent is greater than 1, then always
accept the trial move. If, however, it is less than 1, then accept the
trial move only with a certain probability.

IN:  vector representing x-coordinate of all residues;
     vector representing y-coordinate of all residues;
     energy before and after a trial move is made;
     index of residue subjected to a trial move;
     magnitude of trial movement of residue in x- and y-axes;
     temperature.
OUT: updated vector representing x-coordinate of all residues;
     updated vector representing y-coordinate of all residues.
%}
function [xvals, yvals] = metropolis(xvals, yvals, energyBefore, energyAfter, trial_pos, x_moved, y_moved, temperature)
    helper_exponent = exp( -(energyAfter - energyBefore)/temperature );
    if helper_exponent < 1
        if rand > helper_exponent
            xvals(trial_pos) = xvals(trial_pos) - x_moved;
            yvals(trial_pos) = yvals(trial_pos) - y_moved;
        end
    end
end