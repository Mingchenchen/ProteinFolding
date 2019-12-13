%{
Plot initial polypeptide structure.
Let P residues be blue and H residues be red.

IN:  vector representing hydrophobicity of all residues;
     vector representing x-coordinate of all residues;
     vector representing y-coordinate of all residues.
%}
function plotInitialStruct(hydrophobicity, xvals, yvals)
    scatter(xvals(find(~hydrophobicity)), ...
            yvals(find(~hydrophobicity)), ...
            'ob', 'filled');
    scatter(xvals(find(hydrophobicity)), ...
            yvals(find(hydrophobicity)), ...
            'or', 'filled');
    plot(xvals, yvals, 'color', [0 0.447 0.741]);
    text(xvals(1), yvals(1), ...
         'N', 'HorizontalAlignment', 'right', 'FontSize', 14);
    text(xvals(length(hydrophobicity)), ...
         yvals(length(hydrophobicity)), ...
         'C', 'FontSize', 14);
    xlim([0 length(hydrophobicity)]);
    ylim([-length(hydrophobicity)/2 length(hydrophobicity)/2]);
    axis equal;
    legend('polar residue', 'hydrophobic residue');
    title('Initial Protein State');
    xlabel('x'); ylabel('y');
    set(gca, 'FontSize', 15);
    box on;
end