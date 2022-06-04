N = 1;
folder = "high_accuracy_slack";

for i = 1:N
    load("workspaces/thesis/final_problem3/"+folder+"/iteration_"+num2str(i)+".mat");
    
    if method.save_every_iteration || ~isempty(method.plot_iterations)
        method.plot_name = method.og_plot_name + num2str(iterCount);
    end
    
    displayErrors(errs, M, method.err_treshold, method.err_priority_treshold, method.save_plots || method.save_every_iteration || ismember(i, method.plot_iterations), method.plot_name, method);
end