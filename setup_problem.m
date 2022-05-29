function problem = setup_problem(problem_switch)
    switch problem_switch
        case 0 % chicane
            L1 = 100; L2 = 120; L3 = 160; L4 = 240; % original
            %L1 = 600; L2 = 620; L3 = 660; L4 = 740; % original
            %L1 = 0; L2 = 20; L3 = 20; L4 = 20; % single left turn
            %L1 = -1.0e-5; L2 = 20; L3 = 60; L4 = 60; % no straights
            %L1 = 100; L2 = 100; L3 = 100; L4 = 100; %straight line --> change disconts as well
            params.L1 = L1; params.L2 = L2; params.L3 = L3; params.L4 = L4; 
            myTrack = chicane(params);
            disconts = [myTrack.s1, myTrack.s2, myTrack.s3];
            reference_name_full = 'reference_chicane_N_150.mat';
            reference_name = "reference_chicane";
        case 1 % sine track
            myTrack = smooth_track(0);
            disconts = [];
            reference_name_full = 'reference_smooth_track_N_150.mat';
            reference_name = "reference_smooth_track";
        case 2 % hairpin
            myTrack = hairpin(struct('L1', 100, 'R', 20, 'L2', 60));
            disconts = [myTrack.s1, myTrack.s2];
            reference_name_full = 'reference_hairpin_N_150.mat';
            reference_name = "reference_hairpin";
        case 3 % generic track
            myTrack = generic_track(struct('rho', @(s) 1./(1.0e-5 + 1.0e-1.*exp(-0.1.*(s-50).^2)), 'total_length', 100));
            disconts = [];
            reference_name_full = 'reference_generic_N_150.mat';
            reference_name = "reference_generic";
        case 4 % smooth hairpin
            myTrack = smooth_hairpin(struct('total_length', 70, 'a', 500, 'm', 10, 's1', 20, 's2', 50));
            %myTrack = smooth_hairpin(struct('total_length', 100, 'a', 500, 'm', 15, 's1', 50, 's2', 80));
            disconts = [myTrack.s1, myTrack.s2];
            reference_name_full = 'reference_smooth_hairpin_N_150.mat';
            reference_name = "reference_smooth_hairpin";
        case 5 % circle
            myTrack = circle(struct('total_length', 2*pi*50, 'R', 50));
            disconts = [];
            reference_name_full = 'reference_circle_N_150.mat';
            reference_name = "reference_circle";
        case 6 % zoomed chicane
            load('chicane_zoom_data_compact.mat');
            L1 = 100; L2 = 120; L3 = 160; L4 = 240; % original
            params.L1 = L1; params.L2 = L2; params.L3 = L3; params.L4 = L4;
            N_first = 23; N_end = 46;
            %N_first = 1; N_end = length(res.t);
            s_start = M.s(N_first); s_end = M.s(N_end);
            myTrack = chicane_zoom(s_start, s_end, params); %23 - 46;
            disconts = [myTrack.complete_track.s1, myTrack.complete_track.s2, myTrack.complete_track.s3] - s_start;
            disconts = disconts(disconts > 0);
            disconts = disconts(disconts < s_end- s_start);
            reference_name_full = 'reference_chicane_N_150.mat';
            reference_name = "reference_chicane";
        case 7 %straight
            params = struct('total_length', 100);
            myTrack = straight(params);
            disconts = [];
            reference_name_full = 'reference_chicane_N_150.mat';
            reference_name = 'reference_chicane';
        case 8 % potential of cutting corners
            L1 = 0; L2 = 20; L3 = 20; L4 = 20; % single left turn
            params.L1 = L1; params.L2 = L2; params.L3 = L3; params.L4 = L4; 
            myTrack = chicane(params);
            reference_name_full = 'reference_chicane_N_150.mat';
            reference_name = 'reference_chicane';
            disconts = [];
        case 9 % smooth controls needed at times
            L1 = 0; L2 = 20; L3 = 20; L4 = 20; % single left turn
            params.L1 = L1; params.L2 = L2; params.L3 = L3; params.L4 = L4; 
            myTrack = chicane(params);
            reference_name_full = 'reference_chicane_N_150.mat';
            reference_name = 'reference_chicane';
            disconts = [];
        case 10 % masking effect
            %L1 = 50; L2 = L1+20; L3 = L2+10; % single left turn
            %params.L1 = L1; params.L2 = L2; params.L3 = L3; 
            %myTrack = singleTurn(params);
            myTrack = hairpin(struct('L1', 100, 'R', 20, 'L2', 50));
            reference_name_full = 'reference_chicane_N_150.mat';
            reference_name = 'reference_chicane';
            disconts = [myTrack.s1, myTrack.s2];
        case 11 % final problem
            myTrack = final_problem();
            disconts = [myTrack.Ls(1:end-1)];
            reference_name_full = 'reference_chicane_N_150.mat';
            reference_name = 'reference_chicane';
        %{
        case 6 % no velocity state
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
            problem.max_accel = 20; problem.min_accel = -5;
            problem.roll_off = @(x) exp(-100.*x.^2);
            problem.max_v = 75;
        case 7 % velocity state, but it is constant
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
            problem.max_accel = 20; problem.min_accel = -5;
            problem.roll_off = @(x) exp(-100.*x.^2);
            problem.max_v = 75;
            %}
    end
    rhs = @(x, u, t) [x(3,:).*sin(x(2,:))./get_s_derivative(myTrack, x, t);
                      x(3,:).*tan(u(2,:))./get_s_derivative(myTrack, x, t) - myTrack.evaluate_angle_derivative(t);
                      u(1,:)./get_s_derivative(myTrack, x, t)];
    %problem = struct('t0', 0, 'nx', 3, 'nu', 2, 'x0', [0;0;50], 'rhs', rhs);
    problem = struct('t0', 0, 'nx', 3, 'nu', 2, 'x0', [0;0;50], 'rhs', rhs);
    problem.myTrack = myTrack;
    problem.disconts = disconts;
    problem.xf = [0, NaN, NaN, NaN];
    problem.tf = myTrack.total_length;
    problem.b = @(s) 3;
    problem.max_accel = 20; problem.min_accel = -5;
    %problem.roll_off = @(x) exp(-100.*x.^2);
    %problem.roll_off = @(x) exp(-10.*x.^2);
    %problem.roll_off = @(x) exp(-10.*x.^2).^4;
    problem.roll_off = @(x) 1./(1+65.*x.^2);
    problem.max_v = 50;
    problem.scale = problem.myTrack.total_length;
    problem.problem_switch = problem_switch;
    problem.reference_name = reference_name;
    problem.reference_name_full = reference_name_full;
    
    if problem_switch == 6
        problem.x0 = res.X{N_first}(:,1);
        problem.xf = res.X{N_end}(:,1);
        problem.N_first = N_first; problem.N_end = N_end;
    elseif problem_switch == 7
        problem.roll_off = @(x) ones(size(x));
        problem.max_v = 150;
        problem.xf = [0; 0; 10];
        problem.x0 = [0; 0; 10];
        problem.rhs = @(x, u, t) [x(3,:).*sin(x(2,:))./get_s_derivative(myTrack, x, t);
                    zeros(size(x(3,:)));
                    u(1,:)./get_s_derivative(myTrack, x, t)];
        problem.nu = 1;
        problem.nx = 3;
    elseif problem_switch == 8
        %problem.b = @(s) 0.5 + sqrt(abs(s-15));
        %problem.b = @(s) 0.2 + 0.3.*abs(s-15);
        problem.b = @(s) 4 - 3.*exp(-1.*(s-15).^2);
    elseif problem_switch == 9
        problem.max_v = 20;
        problem.b = @(s) 4;
    elseif problem_switch == 10
        %problem.max_v = 70;
        %problem.x0 = [0; 0; 10];
        problem.roll_off = @(x) 1./(1+100.*x.^2);
    end
end