function [metric, specifics] = get_quality_metric(res, M, rhs, problem)
% returns a quality metric as a weighted sum of different factors
% taken into account is the objective function (final time), the
% integral of the difference with the geometric reference solution and
% the accuracy of the system dynamics
    
    %load('reference_result_N_100.mat');
    %load('reference_mesh_N_100.mat');
    %res_ref = res_previous; M_ref = M_previous;
    
    % objective function
    obj = res.tf;
    disp("obj  = "+num2str(obj));
    
    % compute geometric difference
    %geom_diff = get_geometric_diff(res, M, res_ref, M_ref, problem);
    geom_diff = 0;
    %disp("geom = "+num2str(geom_diff));
    
    % accuracy of system dynamics
    err = log10(get_error_est(res, M, rhs));
    disp("err  = "+num2str(err));
    
    % feasability of path constraints
    feas = get_feas_metric(res, M, problem);
    disp("feas = "+num2str(feas));
    
    metric = obj + geom_diff + err + feas;
    specifics = [obj; geom_diff; err; feas];
end

function diff = get_geometric_diff(res, M, res_ref, M_ref, problem)
    diff = integral(@(s) norm(evaluate_trajectory(res, M, s) - evaluate_trajectory(res_ref, M_ref, s)), 0, problem.myTrack.total_length, 'arrayValued', true);
end

function feas = get_feas_metric(res, M, problem)
    feas = 0;

    svalues = linspace(-1,1,1000);
    x = evaluate_trajectory(res, M, svalues);
    x1 = x(1,:);
    
    for i = 1:length(svalues)
        if abs(x1(i)) > problem.b
            feas = feas + abs(x1(k)) - problem.b;
        end
    end
end