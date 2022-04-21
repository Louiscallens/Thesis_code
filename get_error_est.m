function [err, rels] = get_error_est(res, M, rhs, method, disconts)
    Nb_inter = length(M.s)-1;
    errors = NaN + zeros(Nb_inter, 1);
    rels = NaN + zeros(Nb_inter, 1);
    
    X = res.X; U = res.U;
    
    for k = 1:Nb_inter
        if M.Nu(k) >= 1
            [errors(k), rels(k)] = get_coll_error_est(M.sc{k}, [X{k}, X{k+1}(:,1)], [U{k}, U{k+1}(:,1)], M.Nu(k), rhs, disconts);
        else
            [errors(k), rels(k)] = get_coll_error_est(M.sc{k}, [X{k}, X{k+1}(:,1)], U{k}, M.Nu(k), rhs, disconts);
        end
    end
    
    displayErrors(rels, M, method.err_treshold, method.err_priority_treshold, method.save_plots, method.plot_name);
    
    err = max(rels);
end