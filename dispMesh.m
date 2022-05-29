function dispMesh(M)
    f = figure; f.Position = [597 427.4 560 137.6];
    scatter(M.s, zeros(size(M.s)), 'k', 'filled');
    xoffset = 10;
    xlim([M.s(1)-xoffset, M.s(end)+xoffset]);
    yticks([]);
    ax = gca; ax.YAxis.Visible = 'off';
    xlabel('$s$', 'interpreter', 'latex');
    %saveas(gca, 'figs/thesis/masking/mesh.eps', 'epsc');
end