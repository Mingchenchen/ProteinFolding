%{
Plot initial polypeptide structure.
Let P residues be blue and H residues be red.

IN:  vector representing hydrophobicity of all residues;
     vector representing x-coordinate of all residues;
     vector representing y-coordinate of all residues;
     boolean value for whether we want to label every single residue in plot.
%}
function plotFinalStruct(hydrophobicity, xvals, yvals, plotIdx)
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
    if plotIdx
        pos_idx = 1:length(hydrophobicity);
        for pos = 1:length(hydrophobicity)
            text(xvals(pos), yvals(pos), num2str(pos_idx(pos)));
        end
    end
    xlim([xvals(round(length(hydrophobicity)/2))-round(length(hydrophobicity)/2) xvals(round(length(hydrophobicity)/2))+round(length(hydrophobicity)/2)]);
    ylim([yvals(round(length(hydrophobicity)/2))-round(length(hydrophobicity)/2) yvals(round(length(hydrophobicity)/2))+round(length(hydrophobicity)/2)]);
    axis equal;
    legend('polar residue', 'hydrophobic residue');
    title('Final Protein State');
    xlabel('x'); ylabel('y');
    set(gca, 'FontSize', 15);
    box on;
end