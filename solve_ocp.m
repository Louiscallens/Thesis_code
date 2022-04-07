function results = solve_ocp(M, problem, M_previous, res_previous)
    opti = casadi.Opti();
    Nb_inter = length(M.s)-1;
    
    % create collocation variables
    X = {};
    for i = 1:Nb_inter
        X{i} = opti.variable(problem.nx, M.Nk(i));
    end
    X{Nb_inter+1} = opti.variable(problem.nx, 1);
    
    U = {};
    for i = 1:Nb_inter
        if M.Nu(i) == 0 || M.Nu(i) == 1
            U{i} = opti.variable(problem.nu, 1);
        else
            U{i} = opti.variable(problem.nu, M.Nk(i));
        end
    end
    if M.Nu(end) ~= 0
        U{end+1} = opti.variable(problem.nu,1);
    end
    
    opti = add_initial_final_constraints(opti, problem, X);
    opti = add_path_constraints(opti, problem, X, U);
    opti = add_coll_constraints(opti, problem, M, X, U);
    opti = add_objective(opti, problem, M, X);
    
    if isempty(res_previous.X{1})
        opti = add_initial_initial_guess(opti, M, X, U, problem);
    else
        opti = add_initial_guess(opti, res_previous, M_previous, M, X, U);
    end
    
    % solve OCP   
    opti.solver('ipopt', struct('expand', true));
    %opti.solver('ipopt', struct('expand', true), struct('tol', 1.0e-15));
    opti.callback(@(i) displayTrajectoryX_intermediate(i, M, opti, X, U, problem));
    sol = opti.solve();
    
    results = construct_result(sol, X, U, M, problem);
end

function opti = add_initial_final_constraints(opti, problem, X)
    switch problem.problem_switch
        case {0, 1, 2}
            opti.subject_to(X{1}(:,1) == problem.x0);
            opti.subject_to(X{end}(1) == problem.xf(1));
        otherwise
            opti.subject_to(X{1}(:,1) == problem.x0);
            opti.subject_to(X{end}(:,1) == problem.xf);
    end
end
function opti = add_path_constraints(opti, problem, X, U)
    switch problem.problem_switch
        case 0
            roll_off = @(x) exp(-100.*x.^2);
            
            max_accel = 20; min_accel = -5;
            u = [U{:}];
            opti.subject_to(min_accel.*roll_off(u(2,:))<=u(1,:));
            opti.subject_to(u(1,:)<= max_accel.*roll_off(u(2,:)));
            opti.subject_to(-pi/4 <= u(2,:)); opti.subject_to(u(2,:)<= pi/4);
            
            x = [X{:}]; x1 = x(1,:); x2 = x(2,:); x3 = x(3,:);
            opti.subject_to(-pi/2 < x2 < pi/2);
            opti.subject_to(0 < x3);
            max_v = 75;
            for i = 1:length(X)-1
                %opti.subject_to(X{i}(3,:) <= max_v.*cos((2.0.*U{i}(2,:))).^100);
                opti.subject_to(X{i}(3,:) <= max_v.*roll_off(U{i}(2,:)));
            end
            
            opti.subject_to(x1 <= problem.b); opti.subject_to(x1 >= -problem.b);
        case {1, 2}
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
        tau = M.sc{i};
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
            opti.subject_to(dot_Pi([X{i}, X{i+1}(:,1)], tau(j)) == problem.rhs(X{i}(:,j), uvalues(:,j), tau(j)));
        end
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
    switch problem.problem_switch
        case 0
            x = [X{:}];
            opti.set_initial(x(1,:), 0);
            opti.set_initial(x(2,:), 0);
            opti.set_initial(x(3,:), 0.01);

            for k = 1:length(M.s)-1
                opti.set_initial(U{k}(1,:), 0);
                for j = 1:size(U{k},2)
                    scurr = M.sc{k}(j);
                    opti.set_initial(U{k}(2,j), atan(1/problem.myTrack.evaluate_radius_curvature(scurr)));
                end
            end
        case 1
            x = [X{:}];
            opti.set_initial(x(1,:), 0);
            opti.set_initial(x(2,:), 0);

            for k = 1:length(M.s)-1
                for j = 1:size(U{k},2)
                    scurr = M.sc{k}(j);
                    opti.set_initial(U{k}(1,j), atan(1/problem.myTrack.evaluate_radius_curvature(scurr)));
                end
            end
        case 2
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
function opti = add_initial_guess(opti, res_previous, M_previous, M, X, U)
    Xinit = new_mesh_evaluateX(res_previous.X, M_previous, M);
    for i = 1:length(Xinit)
        opti.set_initial(X{i}(1:2,:), Xinit{i}(1:2,:));
        X3 = Xinit{i}(3,:);
        X3(X3==0) = 1.0e-5;
        opti.set_initial(X{i}(3,:), X3);
    end
        
    Uinit = new_mesh_evaluateU(res_previous.U, M_previous, M);
    for i = 1:length(Uinit)-1
        opti.set_initial(U{i}, Uinit{i});
    end
end