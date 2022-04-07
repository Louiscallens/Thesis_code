function result = get_s_derivative(myTrack, x, t)
    rho = myTrack.evaluate_radius_curvature(t);
    %result = rho./(rho-x(1,:)).*x(3,:).*cos(x(2,:));
    try v = x(3,:); catch v = 0.01; disp("warning: s-derivative --> catch clause activated"); end
    if rho > 1.0e16
        result = v.*cos(x(2,:));
    else
        result = rho./(rho-x(1,:)).*v.*cos(x(2,:));
    end
end