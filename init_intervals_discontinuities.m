function result = init_intervals_discontinuities(disconts, N, total_length)
    disconts = [0, disconts, total_length];
    lengths = disconts(2:end)-disconts(1:end-1);
    Nb_to_divide = N - length(lengths);
    nbs = 1 + floor(Nb_to_divide.*lengths./total_length);
    
    % if there is a left-over interval, add it to the part where the
    % intervals are largest
    while sum(nbs) < N
        widths = lengths./nbs;
        [~, idx] = max(widths);
        nbs(idx) = nbs(idx) + 1;
    end
    
    result= [];
    for i = 1:length(nbs)
        to_add = linspace(disconts(i), disconts(i+1), nbs(i)+1);
        result = [result, to_add(1:end-1)];
    end
    result = [result, total_length];
end