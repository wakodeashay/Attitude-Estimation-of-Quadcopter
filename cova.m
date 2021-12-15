function p = bound(x,s)
    n1=len(x);
    p[1]=s;
    for i:n1
        p[i+1]=3*sqrt(cov(x[1:i+1]));
    end
end