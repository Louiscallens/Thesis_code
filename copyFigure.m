function copyFigure (fig)
    figure(fig); a1 = gca; f = figure(); copyobj(a1, f);
end