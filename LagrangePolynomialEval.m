function res = LagrangePolynomialEval(X,Y,x)
    % Interpolates exactly through data 
    assert(isvector(X));
    N = numel(X);
    switched = false;
    if size(Y,2) ~= N
        Y = Y';
        assert(size(Y,2) == N);
        switched = true;
    end

    res = 0;

    for j=1:N
        p = 1;
        for i=1:N
            if i~=j
               p = p.*(x-X(i))/(X(j)-X(i));
            end
        end
        if ~switched
            res = res+p.*Y(:,j);
        else
            res = res+p.*Y(j,:);
        end
    end

end


