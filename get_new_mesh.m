function Mnew = get_new_mesh(res, M, problem, method, save_plots, plot_name)
% define a new mesh based on specific rules to decide when to increase the
% order of the polynomial or when to split an interval    
    Nb_inter = length(M.s)-1;
    [~, rels] = get_error_est(res, M, problem.rhs, method, problem.disconts);
    to_split = [];
    to_increase = [];
    priority_increase = [];
    active_perf_const = [];
    for k = 1:Nb_inter
        % check slackness of performance constraints
        slack = [mean(res.Yu{k},2); mean(res.Yx{k}(3:end,:),2)];
        active_perf_const = [active_perf_const, slack < method.slack_performance_treshold];
        if min(slack) > method.slack_performance_treshold
            to_split = [to_split, k];
        end
        
        % check slackness of path constraints
        slack = mean(res.Yx{k}(1:2,:),2);
        if min(slack) < method.slack_path_treshold
            to_increase = [to_increase, k];
        
        % check accuracy of the system dynamics
        elseif rels(k) > method.err_treshold
            to_increase = [to_increase, k];
        end
        
        % check extremely low accuracy intervals
        if rels(k) > method.err_priority_treshold
            priority_increase = [priority_increase, k];
        end
    end
        
    knew = 1;
    Mnew = mesh(1, problem.myTrack.total_length, method.Nmin, method.minUDegree, problem.nu);
    splitted = [];
    increased_u = [];
    increased_x = [];
    for k = 1:Nb_inter
        % treat very innacurate intervals first
        if ismember(k, priority_increase)
            if M.Nk(k) + method.Nstep <= method.Nmax
                [M, Mnew, knew] = increase_nb_coll_pts(M, Mnew, k, knew, method);
                increased_x = [increased_x, k];
            else
                [M, Mnew, knew] = split_interval(M, Mnew, k, knew, method);
                splitted = [splitted, k];
            end
            
        % split intervals where no constraint is active
        elseif ismember(k, to_split)
            if M.Nu(k) < 1 && (k ==1 || ismember(k-1, to_split)) && (k == Nb_inter || ismember(k+1, to_split)) % does not exceed linear controls
            %if M.Nu(k) <= 1 && (k ==1 || ismember(k-1, to_split)) && (k == Nb_inter || ismember(k+1, to_split))
                [M, Mnew, knew] = increase_polynomial_order(M, Mnew, k, knew, method);
                increased_u = [increased_u, k];
            else
                [M, Mnew, knew] = split_interval(M, Mnew, k, knew, method);
                splitted = [splitted, k];
            end
            
        % treat those intervals that should increase their number of
        % collocation points
        elseif ismember(k, to_increase)
            if M.Nk(k) + method.Nstep <= method.Nmax
                [M, Mnew, knew] = increase_nb_coll_pts(M, Mnew, k, knew, method);
                increased_x = [increased_x, k];
            else
                [M, Mnew, knew] = split_interval(M, Mnew, k, knew, method);
                splitted = [splitted, k];
            end
        
        else
            % copy the interval
            [M, Mnew, knew] = copy_interval(M, Mnew, k, knew, method);
            
            % prevent masking effect ! knew is alreay updated, so we use
            % the old value !
            %{
            curr_active = active_perf_const(:,k);
            if ~(curr_active(1) || curr_active(2)) % first control input is not limited
                Mnew.Nu(1,knew-1) = min(2, M.Nu(1,k)+1);
                %disp("incraesing polynomial order of u1");
            end
            if ~(curr_active(3) || curr_active(4)) % second control input is not limited
                % increase polynomial order of second control input
                Mnew.Nu(2,knew-1) = min(2, M.Nu(2,k)+1);
                %disp("incraesing polynomial order of u2");
            end
            %}
        end
    end
    
    Mnew.s = [Mnew.s, M.s(end)];
    Mnew.sc = Mnew.add_collocation_times();
    
    displayMeshUpdate(M, splitted, increased_u, increased_x);
end

function [M, Mnew, knew] = increase_polynomial_order(M, Mnew, k, knew, method)
    Mnew.Nk(knew) = M.Nk(k);
	Mnew.s(knew) = M.s(k);
    Mnew.Nu(:,knew) = M.Nu(:,k)+1;
    knew = knew + 1;
end
function [M, Mnew, knew] = split_interval(M, Mnew, k, knew, method)
    Bk = max(2, ceil(M.Nk(k)/method.Nmin));
    ds = (M.s(k+1)-M.s(k))/Bk;
    %if ds < 2.0
    %    [M, Mnew, knew] = copy_interval(M, Mnew, k, knew, method);
    %    return
    %end
    for j = 0:Bk-1
        Mnew.s(knew+j) = M.s(k) + j*ds;
        Mnew.Nk(knew+j) = method.Nmin;
        Mnew.Nu(:,knew+j) = M.Nu(:,k);
    end
    knew = knew + Bk;
end
function [M, Mnew, knew] = increase_nb_coll_pts(M, Mnew, k, knew, method)
    Mnew.Nk(knew) = M.Nk(k) + method.Nstep;
    Mnew.s(knew) = M.s(k);
    Mnew.Nu(:,knew) = M.Nu(:,k);
    knew = knew + 1;
end
function [M, Mnew, knew] = copy_interval(M, Mnew, k, knew, method)
    Mnew.Nk(knew) = M.Nk(k);
    Mnew.s(knew) = M.s(k);
    Mnew.Nu(:,knew) = M.Nu(:,k);
    knew = knew + 1;
end