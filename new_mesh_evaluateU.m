function Ynew = new_mesh_evaluateU(Y, M, Mnew)
% Evaluate the solution Y which defines the coefficients of the lagrange
% polynomials and is defined on mesh M on the new mesh Mnew
    Ynew = {};
    k_previous = -1;
    
    % iterate over new intervals of the new mesh
    for i = 1:length(Mnew.s)-1
        if Mnew.Nu(i) > 1
            Yvalues = zeros(size(Y{1},1), length(Mnew.sc{i})-1);
        else
            Yvalues = zeros(size(Y{1},1), 1);
        end
        
        % iterate over all collocation points within new intervals
        for j = 1:length(Mnew.sc{i})-1
            k = get_interval_idx(M,Mnew.sc{i}(j), k_previous);
            k_previous = k;
            scoll = M.sc{k};
            if Mnew.Nu(i) <= 1
                Yvalues(:,j) = Y{k}(:,1);
                break;
            end
            if M.Nu(k) <= 0
                Yvalues(:,j) = Y{k}(:,1);
            else
                Ycoll = [Y{k}, Y{k+1}(:,1)];
                if size(Ycoll,2) < size(scoll,2)
                    Ycoll = Ycoll(:,1) + [(scoll-scoll(1)).*(Ycoll(1,end)-Ycoll(1,1))./(scoll(end)-scoll(1));(scoll-scoll(1)).*(Ycoll(2,end)-Ycoll(2,1))./(scoll(end)-scoll(1))];
                end
                Yvalues(:,j) = LagrangePolynomialEval(scoll, Ycoll, Mnew.sc{i}(j));
            end
        end
        
        Ynew{i} = Yvalues;
    end
    
    % Store value end tf
    if size(Y{end},2) == 1
        Ynew{end+1} = Y{end}(:,1);
    else
        Ynew{end+1} = LagrangePolynomialEval(M.sc{end}, [Y{end-1}, Y{end}(:,1)], Mnew.s(end));
    end
end