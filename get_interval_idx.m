function k = get_interval_idx(M, s_eval, guess)
    if guess >= 1 && guess <= length(M.s)-1
        if M.s(guess) <= t_eval && t_eval < M.s(guess+1)
            k = guess;
            return
        elseif guess <= length(M.s)-2 && M.s(guess+1) <= s_eval && s_eval < M.s(guess+2)
            k = guess + 1;
            return
        end
    end
    
    k1 = 1; k2 = length(M.s);
    k = k1;
    while k2 > k1+1
        k = k1 + floor((k2-k1)/2);
        if s_eval > M.s(k)
            k1 = k;
        elseif s_eval < M.s(k)
            k2 = k;
        else
            return;
        end        
    end
    k = k1;
end