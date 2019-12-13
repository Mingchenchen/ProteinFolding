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
    
    figure;
    subplot(2, 2, 1); hold on;
    plotInitialStruct(hydrophobicity, xvals, yvals);
    hold off;
    
    for iter = 1:maxiter
        % randomly select a bead/residue
        trial_pos = randi(36);
        
        % randomly select a direction
        % let 1 be N, 2 be S, 3 be E, 4 be W
        direction = randi(4);
        
        % define distance for moving in terms of x- and y-axes
        if direction == 1
            x_moved = 0; y_moved = 1;
        elseif direction == 2
            x_moved = 0; y_moved = -1;
        elseif direction == 3
            x_moved = 1; y_moved = 0;
        elseif direction == 4
            x_moved = -1; y_moved = 0;
        end
        
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
        
        % move selected residue
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
    
    subplot(2, 2, 2); hold on;
    plotFinalStruct(hydrophobicity, xvals, yvals, false);
    hold off;
    
    subplot(2, 2, [3, 4]);
    plotEnergy(energyArray, maxiter);
end