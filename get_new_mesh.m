function [Mnew, updated] = get_new_mesh(res, M, errs, problem, method, save_plots, plot_name)
% define a new mesh based on specific rules to decide when to increase the
% order of the polynomial or when to split an interval    
    updated = false;
    Nb_inter = length(M.s)-1;
    [~, rels] = get_error_est(res, M, problem.rhs, method, problem.disconts);
    to_split = [];
    to_increase = [];
    to_increase_err = [];
    priority_increase = [];
    active_perf_const = [];
    all_slacks = {};
    for k = 1:Nb_inter
        % check slackness of performance constraints
        slack = [mean(res.Yu{k},2); mean(res.Yx{k}(3:end,:),2)];
        all_slacks{end+1} = slack;
        active_perf_const = [active_perf_const, slack < method.slack_performance_treshold];
        if min(slack) > method.slack_performance_treshold
            to_split = [to_split, k];
        end
        
        % check slackness of path constraints
        slack = min(res.Yx{k}(1:2,:),[], 2);
        if min(slack) < method.slack_path_treshold
            to_increase = [to_increase, k];
        
        % check accuracy of the system dynamics
        elseif rels(k) > method.err_treshold
            to_increase_err = [to_increase_err, k];
        end
        
        % check extremely low accuracy intervals
        mean_err = get_weighted_mean(rels, M);
        if rels(k) > method.err_priority_treshold || (rels(k) > method.err_treshold && log10(rels(k)) > mean_err + method.err_order_magnitude_treshold)
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
            %if M.Nk(k) + method.Nstep <= method.Nmax
                [M, Mnew, knew] = increase_nb_coll_pts(M, Mnew, k, knew, method, errs(k));
                increased_x = [increased_x, k];
                updated = true;
            %else
            %    [M, Mnew, knew] = split_interval(M, Mnew, k, knew, 2, method);
            %    splitted = [splitted, k];
            %    updated = true;
            %end
            
        % split intervals where no constraint is active
        elseif ismember(k, to_split) && M.s(k+1)-M.s(k) > method.minimal_interval_width
            if min(M.Nu(:,k)) < 1 && (k ==1 || ismember(k-1, to_split)) && (k == Nb_inter || ismember(k+1, to_split)) % does not exceed linear controls
            %if M.Nu(k) <= 1 && (k ==1 || ismember(k-1, to_split)) && (k == Nb_inter || ismember(k+1, to_split))
                [M, Mnew, knew] = increase_polynomial_order(M, Mnew, k, knew, method);
                increased_u = [increased_u, k];
                updated = true;
            else
                [M, Mnew, knew] = split_interval(M, Mnew, k, knew, 2, method, all_slacks);
                splitted = [splitted, k];
                updated = true;
            end
            
        % treat those intervals that should increase their number of
        % collocation points
        elseif ismember(k, to_increase)
            [M, Mnew, knew] = increase_nb_coll_pts(M, Mnew, k, knew, method);
            increased_x = [increased_x, k];
            updated = true;
        elseif ismember(k, to_increase_err)
            [M, Mnew, knew] = increase_nb_coll_pts(M, Mnew, k, knew, method, errs(k));
            increased_x = [increased_x, k];
            updated = true;
        else
            % copy the interval
            [M, Mnew, knew] = copy_interval(M, Mnew, k, knew, method);
            
            % prevent masking effect ! knew is alreay updated, so we use
            % the old value !
            if method.method_select == 0
                curr_active = active_perf_const(:,k);
                if ~(curr_active(1) || curr_active(2)) % first control input is not limited
                    %Mnew.Nu(1,knew-1) = min(2, M.Nu(1,k)+1);
                    Mnew.Nu(1,knew-1) = max(min(1, M.Nu(1,k)+1), M.Nu(1,k));
                end
                try
                if ~(curr_active(3) || curr_active(4)) % second control input is not limited
                    % increase polynomial order of second control input
                    %Mnew.Nu(2,knew-1) = min(2, M.Nu(2,k)+1);
                    Mnew.Nu(2,knew-1) = max(min(1, M.Nu(2,k)+1), M.Nu(2,k));
                end
                catch
                end
            end
            %}
        end
    end
    
    Mnew.s = [Mnew.s, M.s(end)];
    Mnew.sc = Mnew.add_collocation_times();
    
    displayMeshUpdate(M, splitted, increased_u, increased_x, method);
end

function mean_rels = get_weighted_mean(rels, M)
    % returns average order of magnitude of the errors
    sum = 0;
    for k = 1:length(rels)
        sum = sum + log10(rels(k))*(M.s(k+1)-M.s(k));
    end
    mean_rels = sum/(M.s(end)-M.s(1));
end
function [M, Mnew, knew] = increase_polynomial_order(M, Mnew, k, knew, method)
    Mnew.Nk(knew) = M.Nk(k);
	Mnew.s(knew) = M.s(k);
    Mnew.Nu(:,knew) = min(1,M.Nu(:,k)+1);
    knew = knew + 1;
end
function [M, Mnew, knew] = split_interval(M, Mnew, k, knew, Bk, method, all_slacks)
    if nargin < 7
        ds = (M.s(k+1)-M.s(k))/Bk;
        for j = 0:Bk-1
            Mnew.s(knew+j) = M.s(k) + j*ds;
            Mnew.Nk(knew+j) = max(method.Nmin, ceil(M.Nk(k)/Bk));
            Mnew.Nu(:,knew+j) = M.Nu(:,k);
        end
    else
        ds = get_optimal_cut(all_slacks, M, k, method);
        for j = 0:Bk-1
            Mnew.s(knew+j) = M.s(k) + ds(j+1);
            Mnew.Nk(knew+j) = max(method.Nmin, ceil(M.Nk(k)/Bk));
            Mnew.Nu(:,knew+j) = M.Nu(:,k);
        end
    end
    knew = knew + Bk;
end
function ds = get_optimal_cut(all_slacks, M, k, method)
    curr_slacks = all_slacks{k}(1:2);
    try prev_slacks = all_slacks{k-1}(1:2); catch; prev_slacks = NaN + zeros(size(curr_slacks)); end
    try next_slacks = all_slacks{k+1}(1:2); catch; next_slacks = NaN + zeros(size(curr_slacks)); end
    
    % check if there really is a switch in constraints, otherwise just cut
    % in the middle
    if ~(prev_slacks(1) < method.slack_performance_treshold && next_slacks(2) < method.slack_performance_treshold || ...
            prev_slacks(2) < method.slack_performance_treshold && next_slacks(1) < method.slack_performance_treshold)
         ds = [0, (M.s(k+1)-M.s(k))/2, M.s(k+1)-M.s(k)];
         return
    end
    
    if prev_slacks(1) < curr_slacks(1) % lower bound slack is increasing
        ds = curr_slacks(2)/sum(curr_slacks);
        ds = [ds, ds + curr_slacks(1)/sum(curr_slacks)];
    else                               % upper bound slack is increasing
        ds = curr_slacks(1)/sum(curr_slacks);
        ds = [ds, ds + curr_slacks(2)/sum(curr_slacks)];
    end
    ds = [0, ds];
    ds = ds.*(M.s(k+1)-M.s(k));
end
function [M, Mnew, knew] = increase_nb_coll_pts(M, Mnew, k, knew, method, err)
    if nargin < 6
        NkNew = M.Nk(k) + method.Nstep;
    else
        NkNew = M.Nk(k) + max(method.Nstep, ceil(log(err/method.err_treshold)/log(M.Nk(k))));
    end
    if NkNew > method.Nmax
        Bk = max(2, ceil(NkNew/method.Nmin));
        [M, Mnew, knew] = split_interval(M, Mnew, k, knew, Bk, method);
    else
        %Mnew.Nk(knew) = M.Nk(k) + method.Nstep;
        Mnew.Nk(knew) = NkNew;
        Mnew.s(knew) = M.s(k);
        Mnew.Nu(:,knew) = M.Nu(:,k);
        knew = knew + 1;
    end
end
function [M, Mnew, knew] = copy_interval(M, Mnew, k, knew, method)
    Mnew.Nk(knew) = M.Nk(k);
    Mnew.s(knew) = M.s(k);
    Mnew.Nu(:,knew) = M.Nu(:,k);
    knew = knew + 1;
end