% Newton's method Description
% Starting from a guessed point p0, find the equation of the straight line
% that passes through the point (p0, f(p0)) and has slope of f'(p0).
% The next iterate p1, is the root of this linear equation.

% input:
    % f: input equation
    % p0: starting point
    % tol: tolerance
    % nmax: maximum of iteration
% output: a matrix first row is a vector of iteration; second row is a
%         vector of error per iteration


% The important step is to write a derivative function
function p = newton_method(f, p0, tol, nmax)
    i = 1;
    while i <= nmax
        % write derivative
        % f'(p0) = p0prime
        dx = 0.005;
        p0prime = (feval(f, p0+dx) - feval(f, p0))/dx;
        p = p0 - f(p0)/p0prime;
        
        v(i,1) = i;
        v(i,2) = abs(p-p0);
        v(i,3) = p;
        
        if abs(p-p0) < tol
            p = v;
            break
        else
            p0 = p;
        end
        i = i + 1;
    end
end





