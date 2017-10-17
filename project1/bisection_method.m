% bisection

% input:
    % f: input equation
    % a: min(input interval)
    % b: max(input interval)
    % tol: tolerance
    % nmax: maximum of iteration
% output: a matrix first row is a vector of iteration; second row is a
%         vector of error per iteration

function p = bisection_method(f, a, b, tol, nmax)
    
    if f(a)*f(b) > 0
        disp('f(a) and f(b) have the same sign, there is no root within given interval.');
    else
        i = 1;
        fa = f(a);
        
        while i <= nmax
            p = a + (b-a)/2;
            fp = f(p); % this is err of each iteration
            
            v(i,1) = i; % iteration
            %v(2,i) = a; % a for each step
            %v(3,i) = b; % b for each step
            v(i,2) = abs((b-a)/2); % p for each step, x error
            %v(5,i) = fp; % fp for each step, error
            v(i,3) = p; %root
            
            if (fp == 0) || ((b-a)/2 < tol)
                p = v;
            else
                if fp*fa > 0
                    p = p;
                    a = p;
                    fa = f(p);
                else
                    p = p;
                    b = p;
                end
            end
            i = i + 1;
        end
        %figure
        %plot(vx, vy);
    end
end

