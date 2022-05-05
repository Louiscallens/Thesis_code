function Uex = expandU(U, M)
    Uex = {};
    % loop over the intervals
    for k = 1:length(M.s)-1
        u_toAdd = [];
        % loop over control inputs
        for j = 1:length(U)
            if M.Nu(j,k) == 0
                u_toAdd = [u_toAdd; U{j}{k}.*ones(1, length(M.sc{k}(1:end-1)))];
            elseif M.Nu(j,k) == 1
                u_toAdd = [u_toAdd; U{j}{k} + (M.sc{k}(1:end-1)-M.s(k))./(M.s(k+1)-M.s(k)).*(U{j}{k+1}(1)-U{j}{k}(1))];
            else
                u_toAdd = [u_toAdd; U{j}{k}];
            end
        end
        Uex{k} = u_toAdd;
    end
    
    u_toAdd = [];
    for j = 1:length(U)
        if M.Nu(j,end) ~= 0
            u_toAdd = [u_toAdd; U{j}{end}(end)];
        else
            u_toAdd = [u_toAdd; U{j}{end}];
        end
    end
    %if isempty(u_toAdd(u_toAdd ~= u_toAdd))
        Uex{end+1} = u_toAdd;
    %end
    
end
