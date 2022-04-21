%file1 = 'discont.mat';
%file2 = 'discont_next.mat';

file1 = 'disconts_smooth_hairpin.mat';
file2 = 'disconts_next_smooth_hairpin.mat';

load(file1);
[sfull1, xfull1, dfull1, rfull1, ifull1, scoll1, xcoll1, ucoll1, rcoll1, dcoll1] = generate_plot_values(scoll, xcoll, ucoll, uType, rhs);
load(file2);
[sfull2, xfull2, dfull2, rfull2, ifull2, scoll2, xcoll2, ucoll2, rcoll2, dcoll2] = generate_plot_values(scoll, xcoll, ucoll, uType, rhs);

figure(6); clf;
nx = size(xcoll, 1);
for i = 1:nx
    subplot(nx,4,4*i-3); hold on;
    plot([sfull1, sfull2], [xfull1(i,:), xfull2(i,:)], '-b');
    plot([scoll1, scoll2], [xcoll1(i,:), xcoll2(i,:)], '.b');
    plot([sfull1, sfull2], [ifull1(i,:), ifull2(i,:)], '-r');
    xlabel('$s$', 'interpreter', 'latex');
    ylabel("$x_"+num2str(i)+"$", 'interpreter', 'latex');
    title("accuracy of trajectory", 'interpreter', 'latex');
    xlim([sfull1(1), sfull2(end)]);
    
    subplot(nx,4,4*i-2); hold on;
    plot([sfull1, sfull2], [dfull1(i,:), dfull2(i,:)], '-b');
    plot([scoll1, scoll2], [dcoll1(i,:), dcoll2(i,:)], '.b');
    plot([sfull1, sfull2], [rfull1(i,:), rfull2(i,:)], '-r');
    plot([scoll1, scoll2], [rcoll1(i,:), rcoll2(i,:)], '.r');
    xlabel('$s$', 'interpreter', 'latex');
    ylabel("$x_"+num2str(i)+"$", 'interpreter', 'latex');
    title("comparison of derivative", 'interpreter', 'latex');
    xlim([sfull1(1), sfull2(end)]);
    
    subplot(nx,4,4*i-1);
    semilogy([sfull1, sfull2], abs([xfull1(i,:), xfull2(i,:)] - [ifull1(i,:), ifull2(i,:)]), '-b'); hold on;
    xline([scoll1, scoll2], ':k');
    xlabel('$s$', 'interpreter', 'latex');
    ylabel("$x_"+num2str(i)+"$", 'interpreter', 'latex');
    title("differences in trajectory", 'interpreter', 'latex');
    xlim([sfull1(1), sfull2(end)]);
    
    subplot(nx,4,4*i);
    semilogy([sfull1, sfull2], abs([dfull1(i,:), dfull2(i,:)] - [rfull1(i,:), rfull2(i,:)]), '-b'); hold on;
    xline([scoll1, scoll2], ':k');
    xlabel('$s$', 'interpreter', 'latex');
    ylabel("$x_"+num2str(i)+"$", 'interpreter', 'latex');
    title("differences in derivative", 'interpreter', 'latex');
    xlim([sfull1(1), sfull2(end)]);
end


function [sfull, xfull, dfull, rfull, ifull, scoll, xcoll, ucoll, rcoll, dcoll] = generate_plot_values(scoll, xcoll, ucoll, uType, rhs)
    sfull = linspace(scoll(1), scoll(end), 1000);
    xfull = LagrangePolynomialEval(scoll, xcoll, sfull);
    dfull = LagrangePolynomialDiff(scoll, xcoll, sfull);
    dcoll = LagrangePolynomialDiff(scoll, xcoll, scoll);
    if uType == 0
        ufull = ucoll.*ones(size(sfull));
    elseif uType == 1
        for i = 1:size(ucoll,1)
            ufull(i,:) = ucoll(i,1) + (svalues-svalues(1))/(scoll(end)-scoll(1))*(ucoll(i,end)-ucoll(i,1));
        end
    else
        for i = 1:size(ucoll, 1)
            ufull(i,:) = interp1(scoll(1:end), ucoll(i,:), svalues);
        end
    end
    rfull = rhs(xfull, ufull, sfull);
    rcoll = rhs(xcoll, ucoll, scoll);
    ifull = xcoll(:,1);
    for k = 2:length(sfull)
        ifull(:,k) = ifull(:,k-1) + (rfull(:,k-1)+rfull(:,k))./2.*(sfull(:,k)-sfull(:,k-1));
    end
end