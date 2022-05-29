function [metric, specifics] = get_quality_metric(res, M, rhs, problem, method)
% returns a quality metric as a weighted sum of different factors
% taken into account is the objective function (final time), the
% integral of the difference with the geometric reference solution and
% the accuracy of the system dynamics
    
    %load(problem.reference_name_full);
    %res_ref = res_previous; M_ref = M_previous;
    
    % objective function
    obj = res.tf;
    obj_ref = 0;%res_ref.tf;
    %disp("obj  = "+num2str(obj));
    
    % compute geometric difference
    %geom_diff = get_geometric_diff(res, M, res_ref, M_ref, problem);
    geom_diff = 0;
    %disp("geom = "+num2str(geom_diff));
    
    % accuracy of system dynamics
    err = get_error_est(res, M, rhs, method, problem.disconts);
    %disp("err  = "+num2str(err));
    
    % feasability of path constraints
    feas = get_feas_metric(res, M, problem);
    %disp("feas = "+num2str(feas));

    % number of optimization variables
    nbVars = count_nb_vars(M, problem);
    
    metric = obj + geom_diff + err + feas;
    specifics = [obj; nbVars; err; feas; obj_ref; geom_diff];
end

function diff = get_geometric_diff(res, M, res_ref, M_ref, problem)
    diff = integral(@(s) norm(evaluate_trajectory(res, M, s) - evaluate_trajectory(res_ref, M_ref, s)), 0, problem.myTrack.total_length, 'arrayValued', true);
end

function feas = get_feas_metric(res, M, problem)
    feas = 0;

    %svalues = linspace(-1,1,1000);
    svalues = linspace(0, M.s(end), 1000);
    x = evaluate_trajectory(res, M, svalues);
    x1 = x(1,:);
    
    for i = 1:length(svalues)
        if abs(x1(i)) > problem.b(svalues(i))
            feas = feas + abs(x1(i)) - problem.b(svalues(i));
        end
    end
end

function nbVars = count_nb_vars(M, problem)
    nbVars = 0;
    for k = 1:length(M.s)-1
        % count nb state vars
        nbVars = nbVars + M.Nk(k)*problem.nx;
        
        % count nb of control vars
        for n = 1:problem.nu
            if M.Nu(n,k) == 0 || M.Nu(n,k) == 1
                nbVars = nbVars + 1;
            else
                nbVars = nbVars + M.Nk(k);
            end
        end
    end
end