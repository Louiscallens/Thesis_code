function copyFigureSubplots(fig)
    figure(fig); subplot(211); a1 = gca; subplot(212); a2 = gca; f = figure(); copyobj(a1, f); copyobj(a2,f);
end