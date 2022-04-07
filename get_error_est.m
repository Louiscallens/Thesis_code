function err = get_error_est(res, M, rhs)
    Nb_inter = length(M.s)-1;
    errors = NaN + zeros(Nb_inter, 1);
    rels = NaN + zeros(Nb_inter, 1);
    
    X = res.X; U = res.U;
    
    for k = 1:Nb_inter
        if M.Nu(k) == 1
            [errors(k), rels(k)] = get_coll_error_est(M.sc{k}, [X{k}, X{k+1}(:,1)], [U{k}, U{k+1}], M.Nu(k), rhs);
        else
            [errors(k), rels(k)] = get_coll_error_est(M.sc{k}, [X{k}, X{k+1}(:,1)], U{k}, M.Nu(k), rhs);
        end
    end
    
    err = max(rels);
end