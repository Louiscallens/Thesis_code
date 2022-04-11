function result = init_intervals_discontinuities(disconts, N, total_length)
    disconts = [0, disconts, total_length];
    lengths = disconts(2:end)-disconts(1:end-1);
    nbs = zeros(size(lengths));
    N_left = N;
    for i = 1:length(lengths)-1
        to_add = fix(N_left*lengths(i)/total_length);
        nbs(i) = to_add;
        N_left = N_left - to_add;
    end
    nbs(end) = N_left;
    
    result = [];
    for i = 1:length(nbs)
        to_add = linspace(disconts(i), disconts(i+1), nbs(i)+1);
        result = [result, to_add(1:end-1)];
    end
    result = [result, total_length];
end