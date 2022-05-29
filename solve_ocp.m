function results = solve_ocp(M, problem, problem_switch, method, M_previous, res_previous)
    opti = casadi.Opti();
        
    [X, U, Yx, Yu, U_comp, V] = create_opti_variables(opti, problem, M);
    opti = add_initial_final_constraints(opti, problem, X);
    [opti, Yx, Yu] = add_path_constraints(opti, problem, method, M, X, U, Yx, Yu, V);
    opti = add_coll_constraints(opti, problem, method, M, X, U, V);
    opti = add_objective(opti, method, problem, M, X, U, V);
    
    if isempty(res_previous.X{1})
        opti = add_initial_initial_guess(opti, M, X, U_comp, problem);
    else
        opti = add_initial_guess(opti, res_previous, M_previous, M, X, U_comp, problem, problem_switch);
    end
    
    % solve OCP   
    
    options = struct('nlp_scaling_method','none','mumps_permuting_scaling',0,'mumps_scaling',0,'min_refinement_steps',5,'max_refinement_steps', 100);
    options.max_iter = 10000;
    %options.tol = 1.0e-12;
    %opti.solver('sqpmethod', struct('expand', true, 'convexify_strategy', 'regularize', 'max_iter', 600, ...
    %    'qpsol', 'osqp', 'qpsol_options', struct('osqp', struct('eps_abs', 1.0e-9, 'eps_rel', 1.0e-9))));
    
    %options = struct('nlp_scaling_method','none','linear_solver', 'ma86', 'min_refinement_steps',5,'max_refinement_steps', 100);
    
    opti.solver('ipopt', struct('expand', true), options);
    %opti.callback(@(i) displayTrajectoryX_intermediate(i, M, opti, X, U, Yx, Yu, problem, 50, method));
    try
        sol = opti.solve();
    catch
        sol = opti.debug();
    end
    %figure; spy(sol.value(jacobian(opti.f, opti.x)));
    
    results = construct_result(sol, X, U, Yx, Yu, M, problem);
end

function [X, U, Yx, Yu, U_comp, V] = create_opti_variables(opti, problem, M)
    Nb_inter = length(M.s)-1;
    X = {};
    V = {};
    Yx = {};
    for i = 1:Nb_inter
        X{i} = opti.variable(problem.nx, M.Nk(i));
        V{i} = opti.variable(problem.nx + 2, M.Nk(i));
        %Yx{i} = opti.variable(problem.nx*2, M.Nk(i));
        Yx{i} = casadi.MX.zeros(problem.nx*2, M.Nk(i));
    end
    X{Nb_inter+1} = opti.variable(problem.nx, 1);
    V{Nb_inter+1} = opti.variable(problem.nx + 2, 1);
    %Yx{Nb_inter+1} = opti.variable(problem.nx*2, 1);
    Yx{Nb_inter+1} = casadi.MX.zeros(problem.nx*2, 1);
    
    U = cell(problem.nu, 0);
    Yu = {};
    for i = 1:Nb_inter
        Yu{i} = casadi.MX.zeros(2*problem.nu, M.Nk(i));
        for l = 1:problem.nu
            if M.Nu(l,i) == 0 || M.Nu(l,i) == 1
                U{l}{i} = opti.variable(1, 1);
            else
                U{l}{i} = opti.variable(1, M.Nk(i));
            end
        end
    end
    for l = 1:problem.nu
        if M.Nu(l,end) ~= 0
            U{l}{end+1} = opti.variable(1,1);
        end
    end
    U_comp = U;
    U = expandU(U, M);
    if length(U) > length(Yu)
        Yu{end+1} = opti.variable(problem.nu*2,1);
    end
    
end
function opti = add_initial_final_constraints(opti, problem, X)
    switch problem.problem_switch
        case {0, 1, 2, 3, 4, 5, 8, 9, 10}
            opti.subject_to(X{1}(1:2,1) == problem.x0(1:2));
            opti.subject_to(X{1}(3,1) <= problem.x0(3));
            %opti.subject_to(X{end}(1) == problem.xf(1));
        case 6
            opti.subject_to(X{1}(:,1) == problem.x0);
            %opti.subject_to(X{end} == problem.xf);
        case 7
            opti.subject_to(X{1}(1:2,1) == problem.x0(1:2));
            opti.subject_to(X{end}(1:2,1) == problem.xf(1:2));
            opti.subject_to(X{1}(3,1) <= problem.x0(3));
            opti.subject_to(X{end}(3,1) <= problem.xf(3));
        case 11
            opti.subject_to(X{1}(:,1) == X{end}(:,end));
        otherwise
            opti.subject_to(X{1}(:,1) == problem.x0);
            opti.subject_to(X{end}(:,1) == problem.xf);
    end
end
function [opti, Yx, Yu] = add_path_constraints(opti, problem, method, M, X, U, Yx, Yu, V)
    switch problem.problem_switch
        case {0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11}
            % loop over all intervals
            for k = 1:length(M.s)-1
                get_states = @(s) LagrangePolynomialEval(M.sc{k}, [X{k}, X{k}(:,1)], s);
                get_inputs = @(s) linInterpol(M.sc{k}, [U{k}, U{k+1}(:,1)], s);
                
                % apply constraints to collocation points
                for j = 1:length(M.sc{k})-1
                    %x = get_states(M.sc{k}(j));
                    %u = get_inputs(M.sc{k}(j));
                    x = X{k}(:,j); u = U{k}(:,j);
                    if method.use_viol_vars
                        v = V{k}(:,j);
                        opti.subject_to(problem.min_accel.*problem.roll_off(u(2)) - v(problem.nx + 1).^2 <= u(1) <= problem.max_accel.*problem.roll_off(u(2)) + v(problem.nx + 1).^2);
                        opti.subject_to(2 <= x(3) <= problem.max_v.*problem.roll_off(u(2)) + v(problem.nx + 2).^2);
                        %opti.subject_to(0 <= v);
                    else
                        opti.subject_to(problem.min_accel.*problem.roll_off(u(2)) <= u(1) <= problem.max_accel.*problem.roll_off(u(2)));
                        opti.subject_to(2 <= x(3) <= 2.1+problem.max_v.*problem.roll_off(u(2)));
                    end
                    opti.subject_to(-pi/4 <= u(2) <= pi/4);
                    opti.subject_to(-problem.b(M.sc{k}(j)) <= x(1) <= problem.b(M.sc{k}(j)));
                    opti.subject_to(-pi/2 <= x(2) <= pi/2);
                    
                    Yu{k}(1,j) = U{k}(1,j) - problem.min_accel.*problem.roll_off(U{k}(2,j));
                    Yu{k}(2,j) = problem.max_accel.*problem.roll_off(U{k}(2,j)) - U{k}(1,j);
                    Yu{k}(3,j) = U{k}(2,j) - -pi/4;
                    Yu{k}(4,j) = pi/4 - U{k}(2,j);
                    Yx{k}(1,j) = problem.b(M.sc{k}(j)) - X{k}(1,j);     Yx{k}(2,j) = X{k}(1,j) - -problem.b(M.sc{k}(j));
                    Yx{k}(3,j) = pi/2 - X{k}(2,j);          Yx{k}(4,j) = X{k}(2,j) - -pi/2;
                    Yx{k}(5,j) = X{k}(3,j);                 Yx{k}(6,j) = problem.max_v.*problem.roll_off(U{k}(2,j)) - X{k}(3,j);
                end
                
                % apply constraints to linearly spaced points in the
                % interval
                Nb_extra_constraint_pts = 0;
                svalues = linspace(M.s(k), M.s(k+1), 2+Nb_extra_constraint_pts);
                for j = 2:length(svalues)-1
                    x = get_states(svalues(j));
                    u = get_inputs(svalues(j));
                    opti.subject_to(problem.min_accel.*problem.roll_off(u(2)) <= u(1) <= problem.max_accel.*problem.roll_off(u(2)));
                    opti.subject_to(-pi/4 <= u(2) <= pi/4);
                    opti.subject_to(-problem.b(svalues(j)) <= x(1) <= problem.b(svalues{k}(j)));
                    opti.subject_to(-pi/2 <= x(2) <= pi/2);
                    opti.subject_to(0 <= x(3) <= problem.max_v.*problem.roll_off(u(2)));
                end
            end
            if method.use_viol_vars
                opti.subject_to(problem.min_accel.*problem.roll_off(U{end}(2,1)) - V{end}(problem.nx + 1) <= U{end}(1,1) <= problem.max_accel.*problem.roll_off(U{end}(2,1)) + V{end}(problem.nx + 1));
                opti.subject_to(0 <= X{end}(3,1) <= problem.max_v.*problem.roll_off(U{end}(2,1)) + V{end}(problem.nx + 2));
            else
                opti.subject_to(problem.min_accel.*problem.roll_off(U{end}(2,1)) <= U{end}(1,1) <= problem.max_accel.*problem.roll_off(U{end}(2,1)));
                opti.subject_to(0 <= X{end}(3,1) <= problem.max_v.*problem.roll_off(U{end}(2,1)));
            end
            opti.subject_to(-pi/4 <= U{end}(2,1) <= pi/4);
            opti.subject_to(-problem.b(M.s(end)) <= X{end}(1,1) <= problem.b(M.s(end)));
            opti.subject_to(-pi/2 <= X{end}(2,1) <= pi/2);
            Yu{end}(1,end) = U{end}(1,1) - problem.min_accel.*problem.roll_off(U{end}(2,1));
            Yu{end}(2,end) = problem.max_accel.*problem.roll_off(U{end}(2,1)) - U{end}(1,1);
            Yu{end}(3,end) = U{end}(2,1) - -pi/4;
            Yu{end}(4,end) = pi/4 - U{end}(2,1);
            Yx{end}(1,end) = problem.b(M.s(end)) - X{end}(1,end); Yx{end}(2,end) = X{end}(1,end) - -problem.b(M.s(end));
            Yx{end}(3,end) = pi/2 - X{end}(2,end);      Yx{end}(4,end) = X{end}(2,end) - -pi/2;
            Yx{end}(5,end) = X{end}(3,end);             Yx{end}(6,end) = problem.max_v.*problem.roll_off(U{end}(2,end)) - X{end}(3,end);
            
        case 7
            tol = 1.0e-7;
            
            % loop over all intervals
            for k = 1:length(M.s)-1
                get_states = @(s) LagrangePolynomialEval(M.sc{k}, [X{k}, X{k}(:,1)], s);
                get_inputs = @(s) linInterpol(M.sc{k}, [U{k}, U{k+1}(:,1)], s);
                
                % apply constraints to collocation points
                for j = 1:length(M.sc{k})-1
                    %x = get_states(M.sc{k}(j));
                    %u = get_inputs(M.sc{k}(j));
                    x = X{k}(:,j); u = U{k}(:,j);
                    opti.subject_to(problem.min_accel <= u(1) <= problem.max_accel);
                    %opti.subject_to(-tol <= u(2) <= tol);
                    opti.subject_to(-problem.b(M.sc{k}(j)) <= x(1) <= problem.b(M.sc{k}(j)));
                    opti.subject_to(-pi/2 <= x(2) <= pi/2);
                    
                    opti.subject_to(2 <= x(3) <= problem.max_v);
                    %for t = 1:length(M.sc{k})-1
                    %    if t == j
                    %        uvals = get_inputs(M.sc{k}(t));
                    %        opti.subject_to(0 <= x(3) <= problem.max_v.*problem.roll_off(uvals(2)));
                    %    end
                    %end
                    
                    Yu{k}(1,j) = U{k}(1,j) - problem.min_accel;
                    Yu{k}(2,j) = problem.max_accel - U{k}(1,j);
                    Yx{k}(1,j) = problem.b(M.sc{k}(j)) - X{k}(1,j);     Yx{k}(2,j) = X{k}(1,j) - -problem.b(M.sc{k}(j));
                    Yx{k}(3,j) = pi/2 - X{k}(2,j);          Yx{k}(4,j) = X{k}(2,j) - -pi/2;
                    Yx{k}(5,j) = X{k}(3,j);                 Yx{k}(6,j) = problem.max_v - X{k}(3,j);
                end
                
                % apply constraints to linearly spaced points in the
                % interval
                Nb_extra_constraint_pts = 0;
                svalues = linspace(M.s(k), M.s(k+1), 2+Nb_extra_constraint_pts);
                for j = 2:length(svalues)-1
                    x = get_states(svalues(j));
                    u = get_inputs(svalues(j));
                    opti.subject_to(problem.min_accel <= u(1) <= problem.max_accel);
                    opti.subject_to(-problem.b(svalues(j)) <= x(1) <= problem.b(svalues(j)));
                    opti.subject_to(-pi/2 <= x(2) <= pi/2);
                    opti.subject_to(0 <= x(3) <= problem.max_v);
                end
            end
            opti.subject_to(problem.min_accel <= U{end}(1,1) <= problem.max_accel);
            opti.subject_to(-problem.b(M.s(end)) <= X{end}(1,1) <= problem.b(M.s(end)));
            opti.subject_to(-pi/2 <= X{end}(2,1) <= pi/2);
            opti.subject_to(0 <= X{end}(3,1) <= problem.max_v);
            Yu{end}(1,end) = U{end}(1,1) - problem.min_accel;
            Yu{end}(2,end) = problem.max_accel - U{end}(1,1);
            Yx{end}(1,end) = problem.b(M.s(end)) - X{end}(1,end); Yx{end}(2,end) = X{end}(1,end) - -problem.b(M.s(end));
            Yx{end}(3,end) = pi/2 - X{end}(2,end);      Yx{end}(4,end) = X{end}(2,end) - -pi/2;
            Yx{end}(5,end) = X{end}(3,end);             Yx{end}(6,end) = problem.max_v - X{end}(3,end);
            
            
        otherwise
            u = [U{:}];            
            opti.subject_to(-pi/4 <= u(1,:));
            opti.subject_to(u(1,:)<= pi/4);
            
            x = [X{:}]; x1 = x(1,:); x2 = x(2,:);
            opti.subject_to(-pi/2 < x2 < pi/2);
            opti.subject_to(x1 <= problem.b); opti.subject_to(x1 >= -problem.b);
    end
end
function opti = add_coll_constraints(opti, problem, method, M, X, U, V)
    s = casadi.MX.sym('s');
    for i = 1:length(X)-1
        % collocation
        Xcoll = casadi.MX.sym('Xcoll', problem.nx, M.Nk(i)+1);
        tau = M.sc{i}; %tau_full = [tau, M.s(i+1)];
        polEval  = casadi.Function('polEval', {Xcoll, s}, {LagrangePolynomialEval(tau, Xcoll, s)});
        dot_Pi = casadi.Function('dot_Pi', {Xcoll, s}, {jacobian(polEval(Xcoll, s), s)},{'X', 's'},{'dPi'});
        
        uvalues = [U{i}, U{i+1}(:,1)];
        xvalues = [X{i}, X{i+1}(:,1)];
        
        % apply collocation constraints to the collocation points
        for j = 1:M.Nk(i)
            %if method.use_viol_vars
            %    opti.subject_to(problem.rhs(X{i}(:,j), uvalues(:,j), tau(j)) - V{i}(1:problem.nx,j) <= dot_Pi(xvalues, tau(j)) <= problem.rhs(X{i}(:,j), uvalues(:,j), tau(j)) + V{i}(1:problem.nx,j)); % forward collocation
            %else
                opti.subject_to(dot_Pi(xvalues, tau(j)) == problem.rhs(X{i}(:,j), uvalues(:,j), tau(j))); % forward collocation
            %end
            
            %opti.subject_to(LagrangePolynomialDiff2(tau, xvalues, tau(j)) == problem.rhs(X{i}(:,j), uvalues(:,j), tau(j))); % forward collocation
            %opti.subject_to(dot_Pi(xvalues, tau(j+1)) == problem.rhs(xvalues(:,j+1), uvalues(:,j+1), tau(j+1))); % backward collocation
        end
    end
    % apply collocation constraint to final point as well
	%opti.subject_to(dot_Pi(xvalues, tau(end)) == problem.rhs(X{end}(:,1), uvalues(:,end), tau(end)));
end
function opti = add_objective(opti, method, problem, M, X, U, V)
    switch problem.problem_switch
        otherwise
            int_approx = 0;
            regularization = 0;
            viol_cost = 0;
            %counter = 0;
            for i = 1:length(M.s)-1
                % add time information
                I = get_integration_matrix(M.sc{i});
                int_vals = NaN + casadi.MX.zeros(size(M.sc{i}));
                xvals = [X{i}, X{i+1}(:,1)];
                for j = 1:length(M.sc{i})
                    int_vals(j) = 1/get_s_derivative(problem.myTrack, xvals(:,j), M.sc{i}(j));
                end
                int_approx = int_approx + transpose(I(end,:)*int_vals');
                
                % add regulirization
                %regularization = regularization + sumsqr(U{i} - mean(U{i},2)); % penalize difference from mean in this control interval
                for n = 1:problem.nu
                    if M.Nu(n,i) ~= 0
                        reg = sumsqr(...                % penalize difference from linear control in this interval
                                (U{i}(n,2:end) - (U{i}(n,1) + (M.sc{i}(2:end-1)-M.s(i))./(M.s(i+1)-M.s(i)).*(U{i+1}(n,1)-U{i}(n,1)))));
                        reg = reg/(M.s(i+1)-M.s(i));    % scale inversly with interval width
                        reg = reg/(1.0e1 + (U{i}(n,1)-U{i+1}(n,1))^2);
                        regularization = regularization + reg;
                    end
                end
                %                                ./ (1.0e-5 + sumsqr([U{i}(n,:), U{i+1}(n,1)]))...
                %regularization = regularization + 1/(problem.max_accel-problem.min_accel)*sumsqr((U{i}(1,2:end) - (U{i}(1,1) + (M.sc{i}(2:end-1)-M.s(i))./(M.s(i+1)-M.s(i)).*(U{i+1}(1,1)-U{i}(1,1))))); % penalize difference from linear control in this interval
                %regularization = regularization + 2/pi*sumsqr((U{i}(2,2:end) - (U{i}(2,1) + (M.sc{i}(2:end-1)-M.s(i))./(M.s(i+1)-M.s(i)).*(U{i+1}(2,1)-U{i}(2,1))))); % penalize difference from linear control in this interval
                
                % add constraint violation cost
                %viol_cost = viol_cost + sum(sum(V{i}));
                viol_cost = viol_cost + sumsqr(V{i});
                
                %}
                %{
                for j = 1:length(M.sc{i})-2
                    a = 1/get_s_derivative(problem.myTrack, X{i}(:,j), M.sc{i}(j));
                    b = 1/get_s_derivative(problem.myTrack, X{i}(:,j+1), M.sc{i}(j+1));
                    width = M.sc{i}(j+1)-M.sc{i}(j);
                    int_approx = int_approx + (a+b)/2*width;
                    %counter = counter + 1;
                end
                a = 1/get_s_derivative(problem.myTrack, X{i}(:,end), M.sc{i}(end-1));
                b = 1/get_s_derivative(problem.myTrack, X{i+1}(:,1), M.sc{i}(end));
                width = M.sc{i}(end)-M.sc{i}(end-1);
                int_approx = int_approx + (a+b)/2*width;
                %counter = counter + 1;
                %}
            end
            %disp(counter);
            opti.minimize(int_approx ...
                + method.regularization_weight*regularization ...
                + method.viol_cost_weight*viol_cost);

            %}
            %x = [X{:}];
            %u = [U{:}];
            %opti.minimize(sumsqr(x(1,:)) + sumsqr(u(1,:)));
            %opti.minimize(sumsqr(x(1,:)) - 1.0e-1*sumsqr(x(3,:)));
    end
end
function opti = add_initial_initial_guess(opti, M, X, U_comp, problem)
    switch problem.problem_switch
        case {0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11}
            initial_velocity = 2;
            x = [X{:}];
            opti.set_initial(x(1,:), 0);
            opti.set_initial(x(2,:), 0);
            opti.set_initial(x(3,:), initial_velocity);

            for k = 1:length(M.s)-1
                uvals = atan(1./problem.myTrack.evaluate_radius_curvature(M.sc{k}(1:end-1)));
                
                opti.set_initial(U_comp{1}{k}, 0);
                for j = 1:length(U_comp{2}{k})
                    scurr = M.sc{k}(j);
                    uval = atan(1/problem.myTrack.evaluate_radius_curvature(scurr));
                    opti.set_initial(U_comp{2}{k}(j), uval);
                end
            end
            
        case 7
            initial_velocity = 2;
            x = [X{:}];
            opti.set_initial(x(1,:), 0);
            opti.set_initial(x(2,:), 0);
            opti.set_initial(x(3,:), initial_velocity);

            for k = 1:length(M.s)-1
                uvals = atan(1./problem.myTrack.evaluate_radius_curvature(M.sc{k}(1:end-1)));
                
                opti.set_initial(U_comp{1}{k}, 0);
            end
        case 8
            x = [X{:}];
            opti.set_initial(x(1,:), 0);
            opti.set_initial(x(2,:), 0);
            opti.set_initial(x(3,:), 0.01);

            for k = 1:length(M.s)-1
                for j = 1:size(U{k},2)
                    scurr = M.sc{k}(j);
                    opti.set_initial(U{k}(1,j), atan(1/problem.myTrack.evaluate_radius_curvature(scurr)));
                end
            end
    end 
end
function opti = add_initial_guess(opti, res_previous, M_previous, M, X, U_comp, problem, problem_switch)
    Xinit = new_mesh_evaluateX(res_previous.X, M_previous, M);
    for i = 1:length(Xinit)
        opti.set_initial(X{i}, Xinit{i});
    end
        
    Uinit = new_mesh_evaluateU(res_previous.U, M_previous, M, problem_switch);
    % for case 6 (zoomed chicane), this will fail because the structure of
    % U has changed
    for n = 1:length(U_comp)
        for i = 1:length(Uinit{n})-1
            opti.set_initial(U_comp{n}{i}, Uinit{n}{i});
        end
    end
    
    %[Yxinit, Yuinit] = initialize_slack_variables(M, Xinit, Uinit, problem);
    %for k = 1:length(Yx)
    %    opti.set_initial(Yx{k}, Yxinit{k});
    %end
    %for k = 1:length(Yu)
    %    opti.set_initial(Yu{k}, Yuinit{k});
    %end
    
    %displayTrajectoryX(res_previous, M_previous, problem, false, "");
    %[tc, t] = add_times_to_result(Xinit, M, problem);
    %res = struct('X', {Xinit}, 'U', {Uinit}, 'Yx', {Yxinit}, 'Yu', {Yuinit}, ...
    %    'tc', {tc}, 't', t, 'tf', t(end));
    %displayTrajectoryX(res, M, problem, false, "");
end