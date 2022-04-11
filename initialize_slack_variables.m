function [Yx, Yu] = initialize_slack_variables(M, X, U, problem)
    % add control slackness
    for k = 1:length(U)
        Yu{k}(1,:) = U{k}(1,:) - problem.min_accel.*problem.roll_off(U{k}(2,:));
        Yu{k}(2,:) = problem.max_accel.*problem.roll_off(U{k}(2,:)) - U{k}(1,:);
        Yu{k}(3,:) = U{k}(2,:) - -pi/4;
        Yu{k}(4,:) = pi/4 - U{k}(2,:);
    end
    
    % add state slackness
    for k = 1:length(X)-1
        Yx{k}(1,:) = problem.b - X{k}(1,:);
        Yx{k}(2,:) = X{k}(1,:) - -problem.b;
        Yx{k}(3,:) = pi/2 - X{k}(2,:);
        Yx{k}(4,:) = X{k}(2,:)- -pi/2;
        Yx{k}(5,:) = X{k}(3,:);
        if M.Nu(k) == 0
            uvalues = U{k}(2,:).*ones(size(M.sc{k}(1:end-1)));
        elseif M.Nu(k) == 1
            uvalues = U{k}(2,1) + (M.sc{k}(1:end-1)-M.sc{k}(1))./(M.sc{k}(end)-M.sc{k}(1)).*(U{k+1}(2,1)-U{k}(2,1));
        else
            uvalues = U{k}(2,:);
        end
        Yx{k}(6,:) = problem.max_v.*problem.roll_off(uvalues) - X{k}(3,:);
    end
    
    % add final state slackness
    Yx{end+1}(1,:) = problem.b - X{end}(1,:);
    Yx{end}(2,:)   = X{end}(1,:) - -problem.b;
    Yx{end}(3,:)   = pi/2 - X{end}(2,:);
    Yx{end}(4,:)   = X{end}(2,:)- -pi/2;
    Yx{end}(5,:)   = X{end}(3,:);
    Yx{end}(6,:)   = problem.max_v.*problem.roll_off(U{end}(2,1)) - X{end}(3,1);
end