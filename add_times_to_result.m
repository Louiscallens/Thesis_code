function [tc, t] = add_times_to_result(xres, M, problem)
    tc{1} = convert_s_to_t(M.sc{1}, [xres{1}, xres{2}(:,1)], problem.myTrack, 0);
    t = tc{1}(1);
    for k = 2:length(M.s)-1
        tc{k} = convert_s_to_t(M.sc{k}, [xres{k}, xres{k+1}(:,1)], problem.myTrack, tc{k-1}(end));
        t = [t, tc{k}(1)];
    end
    t = [t, tc{end}(end)];
end