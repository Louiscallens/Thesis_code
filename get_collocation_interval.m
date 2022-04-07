function [tau, w] = get_collocation_interval(t0, tf, d)
    [tau, w] = lgrnodes(d-1);
    tau = tau'; w = w';
    tau = (tf-t0)/2.*(tau-tau(1)) + t0;
    tau = [tau, tf];
end