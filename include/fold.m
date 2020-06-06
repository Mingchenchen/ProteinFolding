%{
Run protein folding.

IN:  vector representing hydrophobicity of all residues;
     vector representing x-coordinate of all residues;
     vector representing y-coordinate of all residues;
     vector representing energy at all iterations of the Metropolis algorithm;
     temperature;
     total number of iterations of Metropolis algorithm.
OUT: updated vector representing x-coordinate of all residues;
     updated vector representing y-coordinate of all residues;
     updated vector represting energy at all iterations of Metropolis algorithm.
%}
function [xvals, yvals, energyArray] = fold(hydrophobicity, xvals, yvals, energyArray, temperature, maxiter)
    energyBefore = 0;  % energy before a given trial move
    energyAfter = 0;   % energy after a given trial move
    
    for iter = 1:maxiter
        % randomly select a bead/residue
        trial_pos = randi(length(hydrophobicity));
        
        % calculate energyBefore for whole polypeptide
        % simplification: consider only part that changes, i.e. trial_pos
        if iter == 1  % else energyAfter will already have been saved as energyBefore
            for pos = 1:length(hydrophobicity)
                if pos ~= trial_pos
                    % energy is sum of Lennard-Jones potential of all pairs of residues
                    potential = calculateEnergyChange(hydrophobicity, xvals, yvals, trial_pos, pos);
                    energyBefore = energyBefore + potential;
                end
            end
            % update energy array for 1st iteration
            energyArray(iter) = energyBefore;
        end
        
        % randomly select a direction in which to move the residue towards
        [x_moved, y_moved] = getDistancesToBeMoved();
        [xvals, yvals] = move(xvals, yvals, trial_pos, x_moved, y_moved);
        
        % calculate energyAfter for whole polypeptide
        % simplification: consider only part that changes, i.e. trial_pos
        for pos = 1:length(hydrophobicity)
            if pos ~= trial_pos
                hasStericClash = checkStericClash(xvals, yvals, trial_pos, pos);
                if hasStericClash
                    [xvals, yvals] = reverseMove(xvals, yvals, trial_pos, x_moved, y_moved);
                    energyAfter = energyBefore;
                    continue;
                    % exits only immediate FOR loop, but the rest of the current iter is essentially O(1).
                    % since energyBefore == energyAfter, reverseMove() will not be triggered again in metropolis().
                end
                
                % energy is sum of Lennard-Jones potential of all pairs of residues
                potential = calculateEnergyChange(hydrophobicity, xvals, yvals, trial_pos, pos);
                energyAfter = energyAfter + potential;
            end
        end
        
        [xvals, yvals] = metropolis(xvals, yvals, energyBefore, energyAfter, trial_pos, x_moved, y_moved, temperature);
        
        % save energyAfter as energyBefore for next iteration
        energyBefore = energyAfter;
        
        % update energy array for all iterations
        energyArray(iter + 1) = energyAfter;
    end
end