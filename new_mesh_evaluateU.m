function Unew = new_mesh_evaluateU(U, M, Mnew)
% Evaluate the solution Y which defines the coefficients of the lagrange
% polynomials and is defined on mesh M on the new mesh Mnew
    Unew = {};
    nu = length(U{1}(:,1));
    k_previous = -1;
    
    % iterate over the control inputs
    for n = 1:nu
        % iterate over new intervals of the new mesh
        for i = 1:length(Mnew.s)-1
            % initialize new control values as zeros
            if Mnew.Nu(n,i) == 0 || Mnew.Nu(n,i) == 1
                k = get_interval_idx(M,Mnew.sc{i}(1), k_previous);
                k_previous = k;
                Uvalues = U{k}(n,1);
            elseif Mnew.Nu(n,i) == 1
                k = get_interval_idx(M,Mnew.sc{i}(1), k_previous);
                k_previous = k;
                Uvalues = linInterpol([M.s(k), M.s(k+1)], [U{k}(n,1), U{k}(n,end)], Mnew.s(i));
            else
                Uvalues = zeros(1, length(Mnew.sc{i})-1);
                % iterate over all points to find the interpolated value
                for j = 1:length(Mnew.sc{i})-1
                    k = get_interval_idx(M,Mnew.sc{i}(j), k_previous);
                    k_previous = k;
                    scoll = M.sc{k};
                    
                    if M.Nu(n,k) == 0
                        Uvalues(j) = U{k}(n,1);
                    elseif M.Nu(n,k) == 1
                        Uvalues(j) = linInterpol([M.s(k), M.s(k+1)], [U{k}(n,1), U{k+1}(n,end)], Mnew.sc{i}(j));
                    else
                        Uvalues(j) = linInterpol(M.sc{k}, [U{k}(n,:), U{k+1}(n,1)], Mnew.sc{i}(j));
                    end
                end
            end
            Unew{n}{i} = Uvalues;
        end
        if size(U{end},2) == 1
            Unew{n}{end+1} = U{end}(n,1);
        else
            Unew{n}{end+1} = linInterpol(M.sc{end}, [U{end-1}(n,:), U{end}(n,1)], Mnew.s(end));
        end
        
    end
end