%{
Imposes two constraints:
 1. two residues cannot share the same pair of x-and y-coordinates,
    because we assume 2D and not 3D.
 2. peptide bond cannot be over- or under-stretched
If either constraint is not met, then the steric clash is resolved by
rejecting the trial move even before the metropolis() function is run.

IN:  vector representing x-coordinate of all residues;
     vector representing y-coordinate of all residues;
     a pair of residues indicated by the indices trial_pos and pos.
OUT: boolean value for whether a steric clash exists.
%}
function hasStericClash = checkStericClash(xvals, yvals, trial_pos, pos)
    hasStericClash = false;
    
    % first constraint
    if xvals(pos) == xvals(trial_pos) && yvals(pos) == yvals(trial_pos)
        hasStericClash = true;
        return;
    end
    
    interResidualDistance = sqrt( (xvals(pos) - xvals(trial_pos)) ^ 2 + ...
                                  (yvals(pos) - yvals(trial_pos)) ^ 2 );
    
    % second constraint
    if pos == trial_pos - 1 || pos == trial_pos + 1
        if interResidualDistance > 1.5 || interResidualDistance < 0.891
            hasStericClash = true;
        end
    end
end