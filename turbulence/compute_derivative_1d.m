% compute the wall-normal derivative of a discretized function, fun(y)
function dfunds = compute_derivative_1d(fun,s)

    dfunds = zeros(size(fun));    % initialize discrete derivative vector

    for j = 1:length(s) % for each index in s
        if j == 1   % if at bottom boundary, use 1-sided derivative
            dfunds(j) = (fun(j+1) - fun(j))/(s(j+1) - s(j));
        elseif j == length(s)% if at top boundary, use 1-sided derivative
            dfunds(j) = (fun(j) - fun(j-1))/(s(j) - s(j-1));
        else    % otherwise, if in the interior of the domain, 2-sided derivative
            dfunds(j) = (fun(j+1) - fun(j-1))/((s(j+1) - s(j-1)));
        end
    end
end