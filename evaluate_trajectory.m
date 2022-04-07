function x = evaluate_trajectory(res, M, s)
    % find the interval we're in
    if length(s) == 1
        idx = get_interval_idx(M, s, 0);
        x = NaN + zeros(size(res.X{1}(:,1),1), length(s));
        for i = 1:size(x,1)
            x(i,:) = LagrangePolynomialEval(M.sc{idx}, [res.X{idx}(i,:), res.X{idx+1}(i,1)], s);
        end
    else
        x = [];
        for k = 1:length(s)
            s_curr = s(k);
            idx = get_interval_idx(M, s_curr, 0);
            for i = 1:size(res.X{1}(:,1),1)
                x(i,k) = LagrangePolynomialEval(M.sc{idx}, [res.X{idx}(i,:), res.X{idx+1}(i,1)], s_curr);
            end
        end
    end
end