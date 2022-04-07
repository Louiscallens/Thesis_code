function displayTrajectoryX_intermediate(iter, M, opti, X, U, problem)
    step = 5;
    if mod(iter, step) ~= 0
        return
    end
    for i = 1:length(X)
        X{i} = opti.debug.value(X{i});
    end
    for i = 1:length(U)
        U{i} = opti.debug.value(U{i});
    end

    curr_res = struct('X', {X}, 'U', {U});
    [tc, t] = add_times_to_result(X, M, problem);
    curr_res.tc = tc; curr_res.t = t; curr_res.tf = t(end);
    
    save_plots = false;
    plot_name = "";
    displayTrajectoryX(curr_res, M, problem, save_plots, plot_name);
    displayTrajectoryU(curr_res, M, problem, save_plots, plot_name);
end