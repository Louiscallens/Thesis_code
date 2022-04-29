function res = LagrangePolynomialDiff2(X,Y,x)
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
        p = 0;
        for i=1:N
            if i~=j
                p_to_add = 1./(X(j)-X(i));
                for m=1:N
                    if m~=j && m~=i
                        p_to_add = p_to_add.*(x-X(m))/(X(j)-X(m));
                    end
                end
                p = p + p_to_add;
            end
        end
        
        if ~switched
            res = res+p.*Y(:,j);
        else
            res = res+p.*Y(j,:);
        end
    end

end


