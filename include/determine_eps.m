%{
Returns the value of \eps for an input pair of residues in a prot.
One residue is located at trial_pos, whilst the other is located at pos.
Protein is input as a one-dimensional vector.

IN: vector representing protein; residue locations in protein indicated by
the indices trial_pos and pos
OUT: value of \eps
%}
function eps = determine_eps(prot, trial_pos, pos)
    % define possible values for eps
    E_HH = -2.3;
    E_HP = -1;
    E_PP = 0;
    E_NN = -10;
    
    % determine value of eps to be used for each pair of residues
    if pos == trial_pos - 1 || pos == trial_pos + 1
        eps = E_NN;
    elseif prot(pos) == 1 && prot(trial_pos) == 1
        eps = E_HH;
    elseif prot(pos) == 0 && prot(trial_pos) == 0
        eps = E_PP;
    elseif prot(pos) ~= prot(trial_pos)
        eps = E_HP;
    end
end