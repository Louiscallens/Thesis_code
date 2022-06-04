function nbVars = count_nb_vars(M, problem)
    nbVars = 0;
    for k = 1:length(M.s)-1
        % count nb state vars
        nbVars = nbVars + M.Nk(k)*problem.nx;
        
        % count nb of control vars
        for n = 1:problem.nu
            if M.Nu(n,k) == 0 || M.Nu(n,k) == 1
                nbVars = nbVars + 1;
            else
                nbVars = nbVars + M.Nk(k);
            end
        end
    end
end