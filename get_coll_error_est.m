function [error_est, rel] = get_coll_error_est(scoll, xcoll, ucoll, uType, rhs)
    [nx, Nk] = size(xcoll);
    nu = size(ucoll,1);
    s0k = scoll(1); sfk = scoll(end);
    
    % get new collocation points
    Mk = Nk + 1;
    s_hat = get_collocation_interval(s0k,sfk,Mk); % gives Mk+1 nodes
    
    % evaluate solution at the new points
    Xs = zeros(nx, Mk+1);
    Us = zeros(nu, Mk+1);
    for i = 1:Mk+1
        for j = 1:nx
            Xs(j,i) = LagrangePolynomialEval(scoll, xcoll(j,:), s_hat(i));
        end
        for j = 1:nu
            if uType == 0
                Us(j,i) = ucoll(j,1);
            elseif uType == 1
                Us(j,i) = ucoll(j,1) + (s_hat(i)-s_hat(1))/(scoll(end)-scoll(1))*(ucoll(j,end)-ucoll(j,1));
            else
                Us(j,i) = LagrangePolynomialEval(scoll(1:end-1), ucoll(j,:), s_hat(i));
            end
        end
    end
    
    % perform integration
    X_hat = zeros(nx, Mk+1);
    I = get_integration_matrix(s_hat);
    
    % SUGGESTION:
    % for i = 2:Mk+1
    %     X_hat(:,i-1) = xcoll(:,1) + transpose(I(i,2:end)*rhs(Xs, Us,
    %     (tf-t0)/2.*t_hat(2:end)+(tf-t0)/2)');
    % end
    
    for i = 1:Mk+1
        X_hat(:,i) = xcoll(:,1) + transpose(I(i,:)*rhs(Xs, Us, s_hat)');
        %X_hat(:,i) = xcoll(:,1) + integral(@(x) full(rhs(LagrangePolynomialEval((tf-t0)/2.*tcoll+(tf+t0)/2, xcoll(j,:), x), LagrangePolynomialEval((tf-t0)/2.*tcoll(1:end-1)+(tf+t0)/2, ucoll(j,1:end-1), x), x)), (tf-t0)/2.*t_hat(1)+(tf+t0)/2, (tf-t0)/2.*t_hat(i)+(tf+t0)/2);
    end    
    
    %error_est = max(abs(X_hat(:,2:end)-Xs(:,2:end)));
    %rel = error_est./(1 + max(Xs,[],2));
    
    error_est = max(max(abs(X_hat(:,2:end)-Xs(:,2:end))));
    rel = max(error_est./(1 + max(Xs,[],2)));
end