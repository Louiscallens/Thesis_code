function y = linInterpol(S, Y, s)
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
    y = Y(:,idx) + (Y(:,idx+1)-Y(:,idx))./(S(idx+1)-S(idx)).*(s-S(idx));
end