% Author: Li-Wei Yap
% Date: 16 October 2018

% to be solved: only one image produced when trying to run multiple protein
% folding sequences concurrently with parfor. Then, I tried adding folding_seq as
% an additional input parameter for the function fold_protein,
% replacing figure with figure(folding_seq), and adding
% saveas(gcf,['prot' num2str(i) '.png']). But sequences end up getting
% printed on the wrong image - perhaps because some sequences finish
% earlier than others? Is there a 'wait' function for parfor?

clear all;  % just a precaution

% re-express the given polypeptide in question as a string of beads
% let 0 be P residue and let 1 be H residue
prot = [0, 0, 1, 0, 1, 1, ...
        0, 1, 0, 0, 0, 0, ...
        1, 1, 0, 0, 1, 1, ...
        0, 0, 0, 1, 0, 1, ...
        1, 0, 0, 0, 1, 1, ...
        0, 0, 1, 1, 0, 1];

% define total no. of iterations of Metropolis algorithm
maxiter = 10^7;  % takes about 150 seconds (2.5 minutes) with 10^7 iterations

% run protein folding
tic
fold_protein(prot, maxiter);
% parfor folding_seq = 1:2
%     fold_protein(prot, maxiter);
% end
toc

function fold_protein(prot, maxiter)
    % define value of temperature
    T = 1;

    % define initial x- and y-values of every residue in polypeptide
    xvals = 1:length(prot);
    yvals = zeros(1, length(prot));
    pos_idx = 1:length(prot);  % use only if you want to label every single residue in plot

    % initialise energy as sum of Lennard-Jones potential of all pairs of residues
    E_i = 0;  % before trial move
    E_f = 0;  % after trial move

    % initialise energy array
    E_array = zeros(1, maxiter + 1);

    % plot initial polypeptide structure
    % let P residues be blue and H residues be red for plot
    % figure(folding_seq);
    figure;
    subplot(2, 2, 1); hold on;
    scatter(xvals(find(~prot)), yvals(find(~prot)), 'ob', 'filled');
    scatter(xvals(find(prot)), yvals(find(prot)), 'or', 'filled');
    plot(xvals, yvals, 'color', [0 0.447 0.741]);
    text(xvals(1), yvals(1), 'N', 'HorizontalAlignment', 'right', 'FontSize', 14);
    text(xvals(length(prot)), yvals(length(prot)), 'C', 'FontSize', 14);
    xlim([0 length(prot)]);
    ylim([-18 18]);
    axis equal;
    legend('polar residue', 'hydrophobic residue');
    title('Initial Protein State');
    xlabel('x'); ylabel('y');
    set(gca, 'FontSize', 15);
    box on; hold off;

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

        % calculate E_i for whole polypeptide
        % simplification: consider only part that changes, i.e. trial_pos
        if iter == 1  % else E_f will already have been saved as E_i (see below)
            % tried parfor for parallelisation, but output not produced after expected amt of time, so terminated early without attempting further
            % maybe it's because prot is only 36 residues long
            for pos = 1:length(prot)
                if pos ~= trial_pos
                    % call helper function (defined below)
                    eps = determine_eps(prot, trial_pos, pos);

                    % determine distance between selected residue and all other residues
                    r = sqrt( (xvals(pos) - xvals(trial_pos)) ^ 2 + ...
                              (yvals(pos) - yvals(trial_pos)) ^ 2 );

                    % calculate Lennard-Jones potential for each pair
                    V = eps * (r^(-12) + 2 * r^(-6));
                    E_i = E_i + V;
                end
            end
            % update energy array for 1st iteration
            E_array(iter) = E_i;
        end

        % move selected residue
        xvals(trial_pos) = xvals(trial_pos) + x_moved;
        yvals(trial_pos) = yvals(trial_pos) + y_moved;

        % calculate E_f for whole polypeptide
        % simplification: consider only part that changes, i.e. trial_pos
        for pos = 1:length(prot)
            if pos ~= trial_pos
                % constraint: two residues cannot share the same pair of x-and y-coordinates, because we assume 2D and not 3D
                if xvals(pos) == xvals(trial_pos) && yvals(pos) == yvals(trial_pos)
                    xvals(trial_pos) = xvals(trial_pos) - x_moved;
                    yvals(trial_pos) = yvals(trial_pos) - y_moved;
                    E_f = E_i;
                    continue;
                end

                % update distance between selected residue and all other residues
                r = sqrt( (xvals(pos) - xvals(trial_pos)) ^ 2 + ...
                          (yvals(pos) - yvals(trial_pos)) ^ 2 );

                % introduce additional constraint for peptide bond (not stated in question)
                if pos == trial_pos - 1 || pos == trial_pos + 1
                    if r > 1.5 || r < 0.891
                        xvals(trial_pos) = xvals(trial_pos) - x_moved;
                        yvals(trial_pos) = yvals(trial_pos) - y_moved;
                        E_f = E_i;
                        continue;
                    end
                end
                
                % call helper function (defined below)
                eps = determine_eps(prot, trial_pos, pos);

                % calculate Lennard-Jones potential for each pair
                V = eps * (r^(-12) + 2 * r^(-6));
                E_f = E_f + V;
            end
        end

        helper_exponent = exp( -(E_f - E_i)/T );
        if helper_exponent < 1  % else accept trial move
            % only accept trial move with a certain probability
            if rand > helper_exponent
                xvals(trial_pos) = xvals(trial_pos) - x_moved;
                yvals(trial_pos) = yvals(trial_pos) - y_moved;
            end
        end

        % save E_f as E_i for next iteration
        % speedup achieved by working with local variables rather than on array itself
        E_i = E_f;

        % update energy array for all iterations
        E_array(iter + 1) = E_f;
    end

    % plot final polypeptide structure
    % let P residues be blue and H residues be red for plot
    subplot(2, 2, 2); hold on;
    scatter(xvals(find(~prot)), yvals(find(~prot)), 'ob', 'filled');
    scatter(xvals(find(prot)), yvals(find(prot)), 'or', 'filled');
    plot(xvals, yvals, 'color', [0 0.447 0.741]);
    text(xvals(1), yvals(1), 'N', 'HorizontalAlignment', 'right', 'FontSize', 14);
    text(xvals(length(prot)), yvals(length(prot)), 'C', 'FontSize', 14);
%     for pos = 1:length(prot)
%         text(xvals(pos), yvals(pos), num2str(pos_idx(pos)));
%     end
    xlim([xvals(18)-18 xvals(18)+18]);
    ylim([yvals(18)-18 yvals(18)+18]);
    axis equal;
    legend('polar residue', 'hydrophobic residue');
    title('Final Protein State');
    xlabel('x'); ylabel('y');
    set(gca, 'FontSize', 15);
    box on; hold off;

    % plot energy array
    subplot(2, 2, [3, 4]);
    plot(E_array);
    set(gca, 'YScale', 'log');
    xlim([0 maxiter]);
    title('Decrease in Energy Over Time Because Protein Folding Is Energetically Favourable');
    xlabel('No. of Iterations'); ylabel('Energy');
    set(gca, 'FontSize', 15);
%     saveas(gcf,['temp' num2str(folding_seq) '.jpg']);
end

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
