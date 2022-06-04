function plot_track(myTrack, b, save_plots, plot_name)
    s_marks = 100:100:myTrack.total_length;
    
    fine_svalues = linspace(0, myTrack.total_length, 1000);
    fine_angles = myTrack.evaluate_angle(fine_svalues);
    [fine_Xs, fine_Ys] = myTrack.evaluate_track_param(fine_svalues);
    plot(fine_Xs, fine_Ys, '--k'); hold on; %center-line
    plot(fine_Xs-b.*sin(fine_angles), fine_Ys+b.*cos(fine_angles), '-k', 'linewidth', 1); % left bound
    plot(fine_Xs+b.*sin(fine_angles), fine_Ys-b.*cos(fine_angles), '-k', 'linewidth', 1); % right bound
    
    axis equal;
    
    angle_marks = myTrack.evaluate_angle(s_marks);
    [Xs_marks, Ys_marks] = myTrack.evaluate_track_param(s_marks);
    for i = 1:length(s_marks)
        plot([Xs_marks(i)+b*sin(angle_marks(i)), Xs_marks(i)-b*sin(angle_marks(i))], ...
             [Ys_marks(i)-b*cos(angle_marks(i)), Ys_marks(i)+b*cos(angle_marks(i))], 'color', [0.8, 0.8, 0.8], 'linewidth', 1);
    end
    plot([0,0], [-b, b], 'k', 'linewidth', 1);
    
    h = annotation('arrow');
    set(h,'parent', gca, ...
        'position', [0,0,8,0], ...
        'HeadLength', 5, 'HeadWidth', 5, 'HeadStyle', 'cback1', ...
        'Color','r','LineWidth',1.5);
    %quiver(0,0,5,0, 'r', 'linewidth', 1.5, 'headlength', 5, 'headwidth', 5);
    
    if save_plots
        %saveas(gca, "figs/case_studies/"+plot_name+"_trajX.eps", 'epsc');
        %saveas(gca, "figs/case_studies/"+plot_name+"_trajX.fig", 'fig');
        %saveas(gca, "figs/case_studies/"+plot_name+"_trajX.png", 'png');
        saveas(gca, plot_name+"_track.eps", 'epsc');
        saveas(gca, plot_name+"_track.fig", 'fig');
        saveas(gca, plot_name+"_track.png", 'png');
    end
end