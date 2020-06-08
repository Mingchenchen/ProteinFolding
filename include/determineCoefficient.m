%{
Returns the value of \eps for an input pair of residues in a prot.

IN:  vector representing hydrophobicity of all residues;
     a pair of residues indicated by the indices trial_pos and pos.
OUT: value of \eps (coefficient in formula).
%}
function coefficient = determineCoefficient(hydrophobicity, trial_pos, pos)
    % define possible values for eps
    E_HH = 2.3;
    E_HP = 1;
    E_PP = 0;
    E_NN = 10;
    
    % determine value of eps to be used for each pair of residues
    if pos == trial_pos - 1 || pos == trial_pos + 1
        coefficient = E_NN;
    elseif hydrophobicity(pos) == 1 && hydrophobicity(trial_pos) == 1
        coefficient = E_HH;
    elseif hydrophobicity(pos) == 0 && hydrophobicity(trial_pos) == 0
        coefficient = E_PP;
    elseif hydrophobicity(pos) ~= hydrophobicity(trial_pos)
        coefficient = E_HP;
    end
end