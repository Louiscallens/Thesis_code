function results = solve_ocp(M, problem, M_previous, res_previous)
    opti = casadi.Opti();
    
    %temp = casadi.MX.sym('temp', 1,1);
    
    [X, U, Yx, Yu] = create_opti_variables(opti, problem, M);
    opti = add_initial_final_constraints(opti, problem, X);
    opti = add_path_constraints(opti, problem, M, X, U, Yx, Yu);
    opti = add_coll_constraints(opti, problem, M, X, U);
    opti = add_objective(opti, problem, M, X);
    
    %[X, U] = create_opti_variables(opti, problem, M);
    %opti = add_initial_final_constraints(opti, problem, X);
    %opti = add_path_constraints(opti, problem, M, X, U);
    %opti = add_coll_constraints(opti, problem, M, X, U);
    %opti = add_objective(opti, problem, M, X);
    
    if isempty(res_previous.X{1})
        %opti = add_initial_initial_guess(opti, M, X, U, Yx, Yu, problem);
        opti = add_initial_initial_guess(opti, M, X, U, problem);
    else
        %opti = add_initial_guess(opti, res_previous, M_previous, M, X, U, Yx, Yu, problem);
        opti = add_initial_guess(opti, res_previous, M_previous, M, X, U, problem);
    end
    
    % solve OCP   
    
    options = struct('nlp_scaling_method','none','mumps_permuting_scaling',0,'mumps_scaling',0,'min_refinement_steps',5,'max_refinement_steps', 20);
    %options = struct('max_iter', 6000);
    %options = struct('tol', 1.0e-15);
    opti.solver('ipopt', struct('expand', true), options);
    %opti.callback(@(i) displayTrajectoryX_intermediate(i, M, opti, X, U, Yx, Yu, problem));
    sol = opti.solve();
    %spy(sol.value(jacobian(opti.f, opti.x)));
    
    results = construct_result(sol, X, U, Yx, Yu, M, problem);
end

function [X, U, Yx, Yu] = create_opti_variables(opti, problem, M)
    Nb_inter = length(M.s)-1;
    X = {};
    Yx = {};
    for i = 1:Nb_inter
        X{i} = opti.variable(problem.nx, M.Nk(i));
        %Yx{i} = opti.variable(problem.nx*2, M.Nk(i));
        Yx{i} = casadi.MX.zeros(problem.nx*2, M.Nk(i));
    end
    X{Nb_inter+1} = opti.variable(problem.nx, 1);
    %Yx{Nb_inter+1} = opti.variable(problem.nx*2, 1);
    Yx{Nb_inter+1} = casadi.MX.zeros(problem.nx*2, 1);
    
    U = {};
    Yu = {};
    for i = 1:Nb_inter
        if M.Nu(i) == 0 || M.Nu(i) == 1
            U{i} = opti.variable(problem.nu, 1);
            %Yu{i} = opti.variable(problem.nu*2, 1);
            Yu{i} = casadi.MX.zeros(problem.nu*2, 1);
        else
            U{i} = opti.variable(problem.nu, M.Nk(i));
            %Yu{i} = opti.variable(problem.nu*2, M.Nk(i));
            Yu{i} = casadi.MX.zeros(problem.nu*2, M.Nk(i));
        end
    end
    if M.Nu(end) ~= 0
        U{end+1} = opti.variable(problem.nu,1);
        Yu{end+1} = opti.variable(problem.nu*2, 1);
    end
end
function opti = add_initial_final_constraints(opti, problem, X)
    switch problem.problem_switch
        case {0, 1, 2, 3, 4, 5}
            opti.subject_to(X{1}(1:2,1) == problem.x0(1:2));
            opti.subject_to(X{1}(3,1) <= problem.x0(3));
            %opti.subject_to(X{end}(1) == problem.xf(1));
        case 6
            opti.subject_to(X{1}(:,1) == problem.x0);
            %opti.subject_to(-1.0e-5 <= X{end} - problem.xf <= 1.0e-5);
            opti.subject_to(X{end} == problem.xf);
        otherwise
            opti.subject_to(X{1}(:,1) == problem.x0);
            opti.subject_to(X{end}(:,1) == problem.xf);
    end
end
function opti = add_path_constraints(opti, problem, M, X, U, Yx, Yu)
    switch problem.problem_switch
        case {0, 1, 2, 3, 4, 5, 6}
            u = [U{:}];
            opti.subject_to(problem.min_accel.*problem.roll_off(u(2,:)) <= u(1,:)<= problem.max_accel.*problem.roll_off(u(2,:)));
            opti.subject_to(-pi/4 <= u(2,:) <= pi/4);
            
            for k = 1:length(U)
                Yu{k}(1,:) =  U{k}(1,:) - problem.min_accel.*problem.roll_off(U{k}(2,:));
                Yu{k}(2,:) = problem.max_accel.*problem.roll_off(U{k}(2,:)) - U{k}(1,:);
                Yu{k}(3,:) = U{k}(2,:) - -pi/4;
                Yu{k}(4,:) = pi/4 - U{k}(2,:);
            end
            
            x = [X{:}]; x1 = x(1,:); x2 = x(2,:); x3 = x(3,:);
            opti.subject_to(-problem.b <= x1 <= problem.b);
            opti.subject_to(-pi/2 <= x2 <= pi/2);
            opti.subject_to(0 <= x3);
            for i = 1:length(X)-1
                if M.Nu(i) ~= 1
                    uvals = U{i}(2,:);
                else
                    uvals = U{i}(2,:) + (M.sc{i}(1:end-1)-M.s(i))./(M.s(i+1)-M.s(i)).*U{i+1}(2,1);
                end
                opti.subject_to(X{i}(3,:) <= problem.max_v.*problem.roll_off(uvals));
                
                Yx{i}(1,:) = problem.b - X{i}(1,:); Yx{i}(2,:) = X{i}(1,:) - -problem.b;
                Yx{i}(3,:) = pi/2 - X{i}(2,:); Yx{i}(4,:) = X{i}(2,:) - -pi/2;
                Yx{i}(5,:) = X{i}(3,:); Yx{i}(6,:) = problem.max_v.*problem.roll_off(uvals) - X{i}(3,:);
            end
            
            %{
            u = [U{:}];
            yu = [Yu{:}];
            opti.subject_to(yu(1,:) == u(1,:) - problem.min_accel.*problem.roll_off(u(2,:))); %opti.subject_to(min_accel.*roll_off(u(2,:))<=u(1,:));
            opti.subject_to(yu(2,:) == problem.max_accel.*problem.roll_off(u(2,:)) - u(1,:)); %opti.subject_to(u(1,:)<= max_accel.*roll_off(u(2,:)));
            opti.subject_to(yu(3,:) == u(2,:) - -pi/4); %opti.subject_to(-pi/4 <= u(2,:));
            opti.subject_to(yu(4,:) == pi/4 - u(2,:));  %opti.subject_to(u(2,:)<= pi/4);
            opti.subject_to(yu(1,:) >= 0);
            opti.subject_to(yu(2,:) >= 0);
            opti.subject_to(yu(3,:) >= 0);
            opti.subject_to(yu(4,:) >= 0);
            
            x = [X{:}]; x1 = x(1,:); x2 = x(2,:); x3 = x(3,:);
            yx = [Yx{:}];
            opti.subject_to(yx(1,:) == problem.b - x1);  %opti.subject_to(x1 <= problem.b);
            opti.subject_to(yx(2,:) == x1 - -problem.b); %opti.subject_to(x1 >= -problem.b);
            opti.subject_to(yx(3,:) == pi/2 - x2);       %opti.subject_to(-pi/2 < x2 < pi/2);
            opti.subject_to(yx(4,:) == x2 - -pi/2);
            opti.subject_to(yx(5,:) == x3);              %opti.subject_to(0 < x3);
            
            for i = 1:length(X)-1
                if M.Nu(i) ~= 1
                    uvals = U{i}(2,:);
                else
                    uvals = U{i}(2,:) + (M.sc{i}(1:end-1)-M.s(i))./(M.s(i+1)-M.s(i)).*U{i+1}(2,1);
                end
                opti.subject_to(Yx{i}(6,:) == problem.max_v.*problem.roll_off(uvals) - X{i}(3,:)); %opti.subject_to(X{i}(3,:) <= max_v.*roll_off(U{i}(2,:)));
            end
            
            opti.subject_to(yx(1,:) >= 0);
            opti.subject_to(yx(2,:) >= 0);
            opti.subject_to(yx(3,:) >= 0);
            opti.subject_to(yx(4,:) >= 0);
            opti.subject_to(yx(5,:) >= 0);
            opti.subject_to(yx(6,:) >= 0);
            %}
            
        otherwise
            u = [U{:}];            
            opti.subject_to(-pi/4 <= u(1,:));
            opti.subject_to(u(1,:)<= pi/4);
            
            x = [X{:}]; x1 = x(1,:); x2 = x(2,:);
            opti.subject_to(-pi/2 < x2 < pi/2);
            opti.subject_to(x1 <= problem.b); opti.subject_to(x1 >= -problem.b);
    end
end
function opti = add_coll_constraints(opti, problem, M, X, U)
    s = casadi.MX.sym('s');
    for i = 1:length(X)-1
        % collocation
        Xcoll = casadi.MX.sym('Xcoll', problem.nx, M.Nk(i)+1);
        tau = M.sc{i}; tau_full = [tau, M.s(i+1)];
        polEval  = casadi.Function('polEval', {Xcoll, s}, {LagrangePolynomialEval(tau, Xcoll, s)});
        dot_Pi = casadi.Function('dot_Pi', {Xcoll, s}, {jacobian(polEval(Xcoll, s), s)},{'X', 's'},{'dPi'});
        
        if M.Nu(i) == 0         % piecewise constant control input
            uvalues = U{i}.*ones(size(U{i},1), length(tau));
        elseif M.Nu(i) == 1     % piecewise linear control input
            uvalues = U{i} + [(tau-tau(1)).*(U{i+1}(1,1)-U{i}(1,1))./(tau(end)-tau(1));(tau-tau(1)).*(U{i+1}(2,1)-U{i}(2,1))./(tau(end)-tau(1))];
        else                    % continuous control input
            uvalues = U{i};
        end
        for j = 1:M.Nk(i)
            opti.subject_to(dot_Pi([X{i}, X{i+1}(:,1)], tau(j)) == problem.rhs(X{i}(:,j), uvalues(:,j), tau(j))); % forward collocation
            %opti.subject_to(dot_Pi([X{i}, X{i+1}(:,1)], tau_full(j+1)) == problem.rhs(X{i}(:,j), uvalues(:,j), tau_full(j+1))); % backward collocation
        end
        % continuity constraint (shouldn't be needed)
        %if i < length(X)-1
        %    opti.subject_to(dot_Pi([X{i}, X{i+1}(:,1)], M.sc{i+1}(1)) == problem.rhs(X{i+1}(:,1), U{i+1}(:,1), M.sc{i+1}(1)));
        %end
    end
end
function opti = add_objective(opti, problem, M, X)
    switch problem.problem_switch
        otherwise
            int_approx = 0;
            for i = 1:length(M.s)-1
                for j = 1:length(M.sc{i})-2
                    a = 1/get_s_derivative(problem.myTrack, X{i}(:,j), M.sc{i}(j));
                    b = 1/get_s_derivative(problem.myTrack, X{i}(:,j+1), M.sc{i}(j+1));
                    width = M.sc{i}(j+1)-M.sc{i}(j);
                    int_approx = int_approx + (a+b)/2*width;
                end
                a = 1/get_s_derivative(problem.myTrack, X{i}(:,end), M.sc{i}(end-1));
                b = 1/get_s_derivative(problem.myTrack, X{i+1}(:,1), M.sc{i}(end));
                width = M.sc{i}(end)-M.sc{i}(end-1);
                int_approx = int_approx + (a+b)/2*width;
            end
            opti.minimize(int_approx);
    end
end
function opti = add_initial_initial_guess(opti, M, X, U, problem)
%function opti = add_initial_initial_guess(opti, M, X, U, Yx, Yu, problem)
    switch problem.problem_switch
        case {0, 1, 2, 3, 4, 5, 6}
            initial_velocity = 1;
            x = [X{:}];
            opti.set_initial(x(1,:), 0);
            opti.set_initial(x(2,:), 0);
            opti.set_initial(x(3,:), initial_velocity);
            
            %yx = [Yx{:}];
            %opti.set_initial(yx(1:2,:), problem.b);
            %opti.set_initial(yx(3:4,:), pi/2);
            %opti.set_initial(yx(5,:), initial_velocity);

            for k = 1:length(M.s)-1
                uvals = atan(1./problem.myTrack.evaluate_radius_curvature(M.sc{k}(1:end-1)));
                %opti.set_initial(Yx{k}(6,:), problem.max_v.*problem.roll_off(uvals) - initial_velocity);
                
                opti.set_initial(U{k}(1,:), 0);
                for j = 1:size(U{k},2)
                    scurr = M.sc{k}(j);
                    uval = atan(1/problem.myTrack.evaluate_radius_curvature(scurr));
                    opti.set_initial(U{k}(2,j), uval);
                    %opti.set_initial(Yu{k}(1,j), 0 - problem.min_accel.*problem.roll_off(uval));
                    %opti.set_initial(Yu{k}(2,j), problem.max_accel.*problem.roll_off(uval) - 0);
                    %opti.set_initial(Yu{k}(3,j), uval - -pi/4);
                    %opti.set_initial(Yu{k}(4,j), pi/4 - uval);
                end
            end
            
        case 7
            x = [X{:}];
            opti.set_initial(x(1,:), 0);
            opti.set_initial(x(2,:), 0);

            for k = 1:length(M.s)-1
                for j = 1:size(U{k},2)
                    scurr = M.sc{k}(j);
                    opti.set_initial(U{k}(1,j), atan(1/problem.myTrack.evaluate_radius_curvature(scurr)));
                end
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
function opti = add_initial_guess(opti, res_previous, M_previous, M, X, U, problem)
    Xinit = new_mesh_evaluateX(res_previous.X, M_previous, M);
    for i = 1:length(Xinit)
        opti.set_initial(X{i}, Xinit{i});
    end
        
    Uinit = new_mesh_evaluateU(res_previous.U, M_previous, M);
    for i = 1:length(Uinit)-1
        opti.set_initial(U{i}, Uinit{i});
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