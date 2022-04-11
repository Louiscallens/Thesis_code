function displayMeshUpdate(M, splitted, increased_u, increased_x)
    split_color = 'b';
    increase_u_color = 'm';
    increase_x_color = 'r';
    nothing_color = 'g';

    f = figure(6); clf; f.Position = [545.0000   71.4000  500.0000  300.0000]; hold on;
    plot(-10,-10, split_color); plot(-10,-10, increase_u_color); 
    plot(-10,-10, increase_x_color); plot(-10,-10, nothing_color);
    legend('$split$', '$increase$ $u$', '$increase$ $nb$ $coll$', '$nothing$',...
        'interpreter', 'latex', 'location', 'best', 'autoupdate', 'off');
    
    for k = 1:length(M.s)-1
        if ismember(k, splitted)
            curr_color = split_color;
        elseif ismember(k, increased_u)
            curr_color = increase_u_color;
        elseif ismember(k, increased_x)
            curr_color = increase_x_color;
        else
            curr_color = nothing_color;
        end
        plot([M.s(k), M.s(k+1)], [0,0], curr_color, 'linewidth', 2);
    end
    xline(M.s, ':k');
    xlim([M.s(1), M.s(end)]);
end