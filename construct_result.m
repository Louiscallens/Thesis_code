function results = construct_result(sol, X, U, Yx, Yu, M, problem)
    xres = {}; ures = {}; yxres = {}; yures = {};
    if M.Nu(end) == 0; stop = length(U); else; stop = length(U)-1; end
    for i = 1:stop
        xres{i} =  full(sol.value(X{i}));  ures{i} =  full(sol.value(U{i}));
        yxres{i} = full(sol.value(Yx{i})); yures{i} = full(sol.value(Yu{i}));
    end
    xres{end+1} = sol.value(X{end}); yxres{end+1} = sol.value(Yx{end});
    ures{end+1} = full(sol.value(U{end})); %yures{end+1} = sol.value(Yu{end});
    
    [tc, t] = add_times_to_result(xres, M, problem);
    
    results = struct('X', {xres}, 'U', {ures}, 'Yx', {yxres}, 'Yu', {yures}, ...
        'tc', {tc}, 't', t, 'tf', t(end));
end