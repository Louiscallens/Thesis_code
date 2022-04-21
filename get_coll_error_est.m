function [error_est, rel] = get_coll_error_est(scoll, xcoll, ucoll, uType, rhs, disconts)
    [nx, Nk] = size(xcoll);
    nu = size(ucoll,1);
    s0k = scoll(1); sfk = scoll(end);
    
    % get new collocation points
    Mk = Nk + 1;
    s_hat = get_collocation_interval(s0k,sfk,Mk); % gives Mk+1 nodes
    
    if ismember(scoll(end), disconts)
        s_hat(end) = s_hat(end) - 1.0e-8;
    end
    
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
                %Us(j,i) = LagrangePolynomialEval(scoll(1:end-1), ucoll(j,:), s_hat(i)); %polynomial interpolation
                Us(j,i) = interp1(scoll(1:end), ucoll(j,:), s_hat(i)); %linear interpolation
            end
        end
    end
    
    % perform integration
    X_hat = zeros(nx, Mk+1);
    I = get_integration_matrix(s_hat);
    
    interpol_rhs = @(s) LagrangePolynomialEval(s_hat, rhs(Xs, Us, s_hat), s);
    
    for i = 1:Mk+1
        X_hat(:,i) = xcoll(:,1) + transpose(I(i,:)*rhs(Xs, Us, s_hat)');
        %X_hat(:,i) = xcoll(:,1) + integral(interpol_rhs, s_hat(1), s_hat(i), 'arrayvalued', true);
        %X_hat(:,i) = xcoll(:,1) + integral(@(s) rhs(LagrangePolynomialEval(scoll, xcoll, s), LagrangePolynomialEval(scoll, ucoll.*ones(size(scoll)), s), s), s_hat(1), s_hat(i), 'arrayvalued', true);
    end    
    
    %error_est = max(abs(X_hat(:,2:end)-Xs(:,2:end)));
    %rel = error_est./(1 + max(Xs,[],2));
    
    error_est = max(abs(X_hat(:,2:end)-Xs(:,2:end)), [], 2);
    rel = max(error_est./(1 + max(Xs,[],2)));
    %rel = max(error_est./max(abs(Xs),[],2));
    error_est = max(error_est);
end