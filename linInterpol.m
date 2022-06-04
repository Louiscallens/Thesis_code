function y = linInterpol(S, Y, svals)
    y = casadi.MX.zeros(size(Y,1), size(svals,2));
    for i = 1:length(svals)
        s = svals(i);
        idx = 1;
        found = false;
        for k = 2:length(S)
            if s < S(k)
                idx = k-1;
                found = true;
                break;
            end
        end
        if ~found
            idx = length(S)-1;
        end
        y(:,i) = Y(:,idx) + (Y(:,idx+1)-Y(:,idx))./(S(idx+1)-S(idx)).*(s-S(idx));
    end
end