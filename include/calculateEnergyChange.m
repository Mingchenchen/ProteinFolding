%{
Returns the value of Lennard-Jones potential for an input pair of residues
in a prot, given that one of them has been subjected to a trial move.

IN:  vector representing hydrophobicity of all residues;
     vector representing x-coordinate of all residues;
     vector representing y-coordinate of all residues;
     a pair of residues indicated by the indices trial_pos and pos.
OUT: value of V (Lennard-Jones potential).
%}
function potential = calculateEnergyChange(hydrophobicity, xvals, yvals, trial_pos, pos)
    coefficient = determineCoefficient(hydrophobicity, trial_pos, pos);  % eps
    
    interResidualDistance = sqrt( (xvals(pos) - xvals(trial_pos)) ^ 2 + ...
                                  (yvals(pos) - yvals(trial_pos)) ^ 2 );  % r
    
    potential = coefficient * (interResidualDistance^(-12) - 2 * interResidualDistance^(-6));  % V
end