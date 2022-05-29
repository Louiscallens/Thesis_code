function [err, rels] = get_error_est(res, M, rhs, method, disconts)
    Nb_inter = length(M.s)-1;
    errors = NaN + zeros(Nb_inter, 1);
    rels = NaN + zeros(Nb_inter, 1);
    
    X = res.X; U = res.U;
    
    for k = 1:Nb_inter
        uvals = [];
        for n = 1:size(U{k},1)
            uvals(n) = (M.Nu(n,k) ~= 1).*U{k}(n,1) + (M.Nu(n,k) == 1).*U{k+1}(n,1);
        end
        uvals = [U{k}, uvals'];
        [errors(k), rels(k)] = get_coll_error_est(M.sc{k}, [X{k}, X{k+1}(:,1)], uvals, rhs, disconts);
    end
    
    displayErrors(rels, M, method.err_treshold, method.err_priority_treshold, method.save_plots, method.plot_name, method);
    
    err = max(rels);
end