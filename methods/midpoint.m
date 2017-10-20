
% f'(x) = (1/2h) * [f(x0+h)-f(x0-h)]

% x = [0: (end-start)/n: end];

% f'(x0) = 1/2h (4 f(x0+h)-3f(x0) - f(x0+2h))


function midpoint = midpoint(f,a,b,n)
    
    h = (b-a)/n;
    xarray = [a:h:b];
    yarray = f(xarray);
    scale = 1/(2*h);
    
    coef = zeros(n+1);
    
    for i = 1: length(xarray)
        
            if i == 1
                
                coef(i,1) = -3*scale;
                coef(i,2) = 4*scale;
                coef(i,3) = -1*scale;
                
            elseif i == length(xarray)
                
                coef(i,end) = -1*scale;
                coef(i,end-1) = 4*scale;
                coef(i,end-2) = -3*scale;
            
            else
                coef(i,i-1) = -1*scale;
                coef(i,i+1) = 1*scale;
    
            end
    end
    
    midpoint = coef * transpose(yarray);
    plot(xarray,mdp,'r')
    hold on
    fplot(f)
    

end





