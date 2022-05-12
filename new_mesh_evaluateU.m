function Unew = new_mesh_evaluateU(U, M, Mnew, problem_switch)
% Evaluate the solution Y which defines the coefficients of the lagrange
% polynomials and is defined on mesh M on the new mesh Mnew
    if problem_switch == 6
        Unew = new_mesh_evaluateU_case6(U, M, Mnew);
        return
    end
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
                        %Uvalues(j) = linInterpol(M.sc{k}(1:end-1), U{k}(n,:), Mnew.sc{i}(j));
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

function Ynew = new_mesh_evaluateU_case6(Y, M, Mnew)
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