function results = construct_result(sol, X, U, M, problem)
    xres = {};
    ures = {};
    if M.Nu(end) == 0; stop = length(U); else; stop = length(U)-1; end
    for i = 1:stop
        xres{i} = full(sol.value(X{i}));
        ures{i} = full(sol.value(U{i}));
    end
    xres{end+1} = sol.value(X{end});
    tau = get_collocation_interval(M.s(end-1), M.s(end), M.Nk(end));
    ures{end+1} = full(sol.value(U{end}));
    
    [tc, t] = add_times_to_result(xres, M, problem);
    
    results = struct('X', {xres}, 'U', {ures}, 'tc', {tc}, 't', t, 'tf', t(end));
end