function times = convert_s_to_t(svalues, x, myTrack, start_time)
    times = NaN + zeros(size(svalues));
    times(1) = start_time;
    for i = 2:length(svalues)
        times(i) = times(i-1) + integral(@(time) 1./get_s_derivative(myTrack, interp1(svalues, x', time)', time), svalues(i-1), svalues(i), 'arrayvalued', 1);
    end
end