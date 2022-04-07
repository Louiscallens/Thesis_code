function Mnew = get_new_mesh(res, M, problem, method, save_plots, plot_name)
% define a new mesh based on specific rules to decide when to increase the
% order of the polynomial or when to split an interval
    Nb_inter = length(M.s)-1;
    
    Mnew = struct('s', [], 'Nk', []);

    knew = 1;
    
    all_cs = NaN + zeros(10, Nb_inter);
    vel_cs = [];
    to_split = [];
    to_increase = [];
    
    figure(5); hold off;
    for k = 1:Nb_inter
        curr_input = res.U{k};
        curr_state = res.X{k};
        mean_state = mean(curr_state, 2);
        
        
        max_accel = 20;
        min_accel = -5;
        max_v = 75;
        
        c = NaN + zeros(10,1);
        %c = NaN + zeros(6,1);
        c(1) = curr_input(1) - min_accel*(1.1*cos(curr_input(2))^4);
        c(2) = max_accel*(1.1*cos(curr_input(2)))^4 - curr_input(1);
        c(3) = curr_input(2) + pi/2;
        c(4) = pi/2 - curr_input(2);
        c(5) = mean_state(3);
        vel_cs = [vel_cs, max_v*cos(2*curr_input(2))^100 - curr_state(3,:)];
        c(6) = mean(max_v*cos(2*curr_input(2))^100 - curr_state(3,:));
        
        for i = 1:length(problem.r)
            c(6+i) = min((curr_state(1,:)-problem.rx(i)).^2+(curr_state(2,:)-problem.ry(i)).^2 - problem.r(i)^2);
        end
        
        all_cs(:,k) = c;
        if min(c(1:6)) < 1.0e-5
            color = 'r';
        elseif min(c(7:end)) < 1.0e-3
            color = 'b';
        else
            color = 'g';
        end
        scatter(k, 0, color); hold on;
        
        % process the interval if no constraint is active
        if min(c(1:6)) > 1.0e-5
            to_split = [to_split, k];
        elseif min(c(7:end)) < 1.0e-3
            to_increase = [to_increase, k];
        end
    end
    for k = 1:Nb_inter
        if ismember(k, to_split)
            % process the interval
            
            % check neighbours
            %if M.Nu(k) < 1 && (k ==1 || ismember(k-1, to_split)) && (k == Nb_inter || ismember(k+1, to_split)) % does not exceed linear controls
            if M.Nu(k) <= 1 && (k ==1 || ismember(k-1, to_split)) && (k == Nb_inter || ismember(k+1, to_split))
                % change to piecewise linear (or higher order)
                Mnew.Nk(knew) = M.Nk(k);
                Mnew.t(knew) = M.t(k);
                Mnew.Nu(knew) = M.Nu(k)+1;
                knew = knew + 1;
            else
                % split the interval
                Bk = 2;
                dt = (M.t(k+1)-M.t(k))/Bk;
                for j = 0:Bk-1
                    Mnew.t(knew+j) = M.t(k) + j*dt;
                    Mnew.Nk(knew+j) = method.Nmin;
                    Mnew.Nu(knew+j) = M.Nu(k);
                end
                knew = knew + Bk;
            end
        elseif ismember(k, to_increase)
            % increase order of approximating polynomial
            %disp("Increasing order");
            Mnew.Nk(knew) = M.Nk(k) + 5;
            Mnew.t(knew) = M.t(k);
            Mnew.Nu(knew) = M.Nu(k);
            knew = knew + 1;
        else % copy the interval to the new mesh
            Mnew.Nk(knew) = M.Nk(k);
            Mnew.t(knew) = M.t(k);
            Mnew.Nu(knew) = M.Nu(k);
            knew = knew + 1;
        end
    end
    
    Mnew.t = [Mnew.t, M.t(end)];
    Mnew.tc = get_collocation_times(Mnew);
    
    if save
        saveas(gca, plot_name+"_mesh_eval.eps", 'epsc');
        saveas(gca, plot_name+"_mesh_eval.fig", 'fig');
        saveas(gca, plot_name+"_mesh_eval.png", 'png');
    end
    
    figure(7);
    subplot(2,1,1); hold off;
    for k = 1:6
        semilogy(all_cs(k,:), 'o-'); hold on;
    end
    title('slackness to performance constraints', 'interpreter', 'latex');
    xlabel('collocation points', 'interpreter', 'latex');
    subplot(2,1,2); hold off;
    for k = 7:size(all_cs,1)
        semilogy(all_cs(k,:), 'o-'); hold on;
    end
    title('slackness to path constraints', 'interpreter', 'latex');
    xlabel('collocation points', 'interpreter', 'latex');
    
    if save
        saveas(gca, plot_name+"_slack.eps", 'epsc');
        saveas(gca, plot_name+"_slack.fig", 'fig');
        saveas(gca, plot_name+"_slack.png", 'png');
    end
    
    figure(10); clf;
    semilogy(vel_cs, 'o-');
    title('velocity slackness at collocation points');
    xlabel('collocation points', 'interpreter', 'latex');
    if save
        saveas(gca, plot_name+"_vel_slack.eps", 'epsc');
        saveas(gca, plot_name+"_vel_slack.fig", 'fig');
        saveas(gca, plot_name+"_vel_slack.png", 'png');
    end
end