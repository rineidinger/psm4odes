function y = serieseval(coefs,t0,t)
%SERIESEVAL Horner's rule to evaluate polynomial(s) at t(s)
% using coefs of increasing powers of (t-t0), t may be a (row) vector,
% each row of coefs is one truncated series (polynomial).
% Output y is column vector(s) for each time in t.
    n = length(coefs);
    val = coefs(:,n);
    h = t-t0;
    for m = n-1:-1:1
        val = val.*h + coefs(:,m);
    end
    y = val;
end