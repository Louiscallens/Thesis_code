function problem = setup_problem(problem_switch)
    switch problem_switch
        case 0 % full problem
            L1 = 100; L2 = 120; L3 = 160; L4 = 240; % original
            %L1 = 0; L2 = 20; L3 = 20; L4 = 20; % single left turn
            %L1 = -1.0e-5; L2 = 20; L3 = 60; L4 = 60; % no straights
            params.L1 = L1; params.L2 = L2; params.L3 = L3; params.L4 = L4; 
            myTrack = track(params);
            rhs = @(x, u, t) [x(3,:).*sin(x(2,:))./get_s_derivative(myTrack, x, t);
                              x(3,:).*tan(u(2,:))./get_s_derivative(myTrack, x, t) - myTrack.evaluate_angle_derivative(t);
                              u(1,:)./get_s_derivative(myTrack, x, t)];
            problem = struct('t0', 0, 'nx', 3, 'nu', 2, 'x0', [0;0;0.01], 'rhs', rhs);
            problem.myTrack = myTrack;
            problem.xf = [0, NaN, NaN, NaN];
            problem.L1 = L1; problem.L2 = L2; problem.L3 = L3; problem.L4 = L4;
            problem.tf = myTrack.total_length;
            problem.b = 4;
        case 1 % no velocity state
            %L1 = 100; L2 = 120; L3 = 160; L4 = 240; % original
            %L1 = -1.0e-5; L2 = 20; L3 = 20; L4 = 20; % single left turn
            %L1 = 0; L2 = 0; L3 = 40; L4 = 40; % single right turn
            L1 = -1.0e-5; L2 = 20; L3 = 60; L4 = 60; % no straights
            params.L1 = L1; params.L2 = L2; params.L3 = L3; params.L4 = L4; 
            myTrack = track(params);
            rhs = @(x, u, t) [0.01.*sin(x(2,:))./get_s_derivative(myTrack, x, t);
                              0.01.*tan(u(1,:))./get_s_derivative(myTrack, x, t) - myTrack.evaluate_angle_derivative(t)];
            problem = struct('t0', 0, 'nx', 2, 'nu', 1, 'x0', [0;0], 'rhs', rhs);
            problem.myTrack = myTrack;
            problem.xf = [0, NaN, NaN, NaN];
            problem.L1 = L1; problem.L2 = L2; problem.L3 = L3; problem.L4 = L4;
            problem.tf = myTrack.total_length;
            problem.b = 4;
        case 2 % velocity state, but it is constant
            %L1 = 100; L2 = 120; L3 = 160; L4 = 240; % original
            L1 = 0; L2 = 20; L3 = 20; L4 = 20; % single left turn
            %L1 = 0; L2 = 0; L3 = 40; L4 = 40; % single right turn
            %L1 = 0; L2 = 20; L3 = 60; L4 = 60; % no straights
            params.L1 = L1; params.L2 = L2; params.L3 = L3; params.L4 = L4; 
            myTrack = track(params);
            rhs = @(x, u, t) [x(3,:).*sin(x(2,:))./get_s_derivative(myTrack, x, t);
                              x(3,:).*tan(u(1,:))./get_s_derivative(myTrack, x, t) - myTrack.evaluate_angle_derivative(t);
                              0];
            problem = struct('t0', 0, 'nx', 3, 'nu', 1, 'x0', [0;0;0.01], 'rhs', rhs);
            problem.myTrack = myTrack;
            problem.xf = [0, NaN, NaN, NaN];
            problem.L1 = L1; problem.L2 = L2; problem.L3 = L3; problem.L4 = L4;
            problem.tf = myTrack.total_length;
            problem.b = 4;
    end
    problem.scale = problem.myTrack.total_length;
    problem.problem_switch = problem_switch;
end