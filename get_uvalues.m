function [tvalues, uvalues, tVars, uVars] = get_uvalues(res, M)
    plot_polynomials = false;
    
    tvalues = [];
    uvalues = [];
    tVars = 0;
    uVars = res.U{1}(:,1);
    for k = 1:length(res.U)-1
        if M.Nu(k) == 0
            tvalues = [tvalues, res.t(k), res.t(k+1)];
            uvalues = [uvalues, res.U{k}, res.U{k}];
            tVars = [tVars, (res.t(k)+res.t(k+1))/2];
            uVars = [uVars, res.U{k}];
        elseif M.Nu(k) == 1
            uvalues = [uvalues, res.U{k} + (res.tc{k}(1:end-1)-res.tc{k}(1)).*(res.U{k+1}(:,1)-res.U{k})./(res.tc{k}(end)-res.tc{k}(1))];
            tvalues = [tvalues, res.tc{k}(1:end-1)];
            tVars = [tVars, res.t(k), res.t(k+1)];
            uVars = [uVars, res.U{k}, res.U{k+1}(:,1)];
        else
            if plot_polynomials
                temp = linspace(res.tc{k}(1), res.tc{k}(end), 100);
                uvalues = [uvalues, LagrangePolynomialEval(res.tc{k}(1:end-1), res.U{k}, temp)];
                tvalues = [tvalues, temp];
                tVars = [tVars, res.tc{k}(1:end-1)];
                uVars = [uVars, res.U{k}];
            else
                uvalues = [uvalues, res.U{k}];
                tvalues = [tvalues, res.tc{k}(1:end-1)];
                tVars = [tVars, res.tc{k}(1:end-1)];
                uVars = [uVars, res.U{k}];
            end
        end
    end
    uvalues = [uvalues, res.U{end}];
    tvalues = [tvalues, res.t(end)];
end