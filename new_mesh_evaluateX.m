function Xnew = new_mesh_evaluateX(X, M, Mnew)
% Evaluate the solution Y which defines the coefficients of the lagrange
% polynomials and is defined on mesh M on the new mesh Mnew
    Xnew = {};
    k_previous = -1;
    
    % iterate over new intervals of the new mesh
    for i = 1:length(Mnew.s)-1
        Xvalues = zeros(size(X{1},1), length(Mnew.sc{i})-1);
        
        % iterate over all collocation points within new intervals
        for j = 1:length(Mnew.sc{i})-1
            k = get_interval_idx(M,Mnew.sc{i}(j), k_previous);
            k_previous = k;
            scoll = M.sc{k};
            Xcoll = [X{k}, X{k+1}(:,1)];
            Xvalues(:,j) = LagrangePolynomialEval(scoll, Xcoll, Mnew.sc{i}(j));
        end
        
        Xnew{i} = Xvalues;
    end
    
    % Store value end tf
    Xnew{end+1} = LagrangePolynomialEval(M.sc{end}, [X{end-1}, X{end}(:,1)], Mnew.s(end));
end