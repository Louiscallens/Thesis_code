function dispMesh(M)
    f = figure; %f.Position = [597 427.4 560 137.6];
    f.Position = [221.8 427.4 935.2 137.6];
    scatter(M.s, zeros(size(M.s)), 15, 'k', 'filled');
    xoffset = 10;
    xlim([M.s(1)-xoffset, M.s(end)+xoffset]);
    yticks([]);
    ax = gca; ax.YAxis.Visible = 'off';
    xlabel('$s$', 'interpreter', 'latex');
    %saveas(gca, 'figs/thesis/masking/mesh.eps', 'epsc');
    saveas(gca, 'figs/thesis/final_problem2/mesh.eps', 'epsc');
end