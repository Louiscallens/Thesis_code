function Ynew = new_mesh_evaluateX(Y, M, Mnew)
% Evaluate the solution Y which defines the coefficients of the lagrange
% polynomials and is defined on mesh M on the new mesh Mnew
    Ynew = {};
    k_previous = -1;
    
    % iterate over new intervals of the new mesh
    for i = 1:length(Mnew.s)-1
        Yvalues = zeros(size(Y{1},1), length(Mnew.tsc{i})-1);
        
        % iterate over all collocation points within new intervals
        for j = 1:length(Mnew.sc{i})-1
            k = get_interval_idx(M,Mnew.sc{i}(j), k_previous);
            k_previous = k;
            scoll = M.sc{k};
            Ycoll = [Y{k}, Y{k+1}(:,1)];
            Yvalues(:,j) = LagrangePolynomialEval(scoll, Ycoll, Mnew.sc{i}(j));
        end
        
        Ynew{i} = Yvalues;
    end
    
    % Store value end tf
    Ynew{end+1} = LagrangePolynomialEval(M.sc{end}, [Y{end-1}, Y{end}(:,1)], Mnew.s(end));
end