%{
Plot energy for all iterations of the Metropolis algorithm in the log-scale.

IN:  vector representing energy at all iterations of the Metropolis algorithm;
     total number of iterations of Metropolis algorithm.
%}
function plotEnergy(energyArray, maxiter)
    plot(energyArray, 'LineWidth', 3);
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xlim([0 maxiter]);
    title('Decrease in Energy Over Time Because Protein Folding Is Energetically Favourable');
    xlabel('No. of Iterations'); ylabel('Energy');
    set(gca, 'FontSize', 15);
end