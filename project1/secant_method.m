% secant method

% input:
    % f - input function that you want to get answer
    % p0 - initial approximates 0
    % p1 - initial approximates 1
    % tol - tolerance
    % nmax - maximum iterations
% output: a matrix first row is a vector of iteration; second row is a
%         vector of error per iteration



function p = secant_method(f, p0, p1, tol, nmax)
    i = 2;
    q0 = f(p0);
    q1 = f(p1);

    while i <= nmax
       p = p1 - q1*(p1-p0)/(q1-q0);
       q = f(p);
       
       v(i-1,1) = i-1;
       v(i-1,2) = abs(p1-p); %err
       v(i-1,3) = p;
       if abs(p1-p) < tol
           p = v;
           break;
       else
           i = i + 1;
           p0 = p1;
           q0 = q1;
           p1 = p;
           q1 = f(p);

       end
    end

end