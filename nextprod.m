% Direct ports of Julia's `nextprod` [1], with helper function `nextpow` [2],
% which are in Julia 1.2.0.
%
% [1] https://github.com/JuliaLang/julia/blob/c6da87ff4bc7a855e217856757ad3413cf6d1f79/base/combinatorics.jl%L248-L262
% [2] https://github.com/JuliaLang/julia/blob/c6da87ff4bc7a855e217856757ad3413cf6d1f79/base/intfuncs.jl%L334-L356
%
% This derivative work is licensed under the MIT License, the same license as Julia Base.
function nprod = nextprod(a, x)
%   """Next integer greater than or equal to `x` that can be written as ``\\prod k_i^{a_i}`` for integers
%   ``a_1``, ``a_2``, etc.
%   % Examples
%   ```jldoctest
%   julia> nextprod([2, 3], 105)
%   108
%   julia> 2^2 * 3^3
%   108
%   ```
%   """
k = length(a);
v = ones(k,1); %[1] * k  % current value of each counter
mx = zeros(k,1);
for q=1:k
    mx(q) = nextpow(a(q), x);  % maximum value of each counter
end
v(1) = mx(1);  % start at first case that is >= x
p = mx(1);  % initial value of product in this case
best = p;
icarry = 2;

while v(end)<mx(end)
    if p >= x
        if p<best
            best=p;
        end
        carrytest = true;
        while carrytest
            p = floor(p / v(icarry - 1));
            v(icarry - 1) = 1;
            icarry = icarry + 1;
            p = p * a(icarry - 1);
            v(icarry - 1) = v(icarry - 1)*a(icarry - 1);
            carrytest = (v(icarry - 1) > mx(icarry - 1)) && icarry < k;
        end
        if p < x
            icarry = 2;
        end
    else
        while p < x
            p = p*a(1);
            v(1) = v(1)*a(1);
        end
    end
end
if mx(end)<best
    nprod=mx(end);
else
    nprod=best;
end
end

function outp = nextpow(a, x)
%   """The smallest `a^n` not less than `x`, where `n` is a non-negative integer.
%
%   `a` must be greater than 1, and `x` must be greater than 0.
%   % Examples
%   ```jldoctest
%   julia> nextpow(2, 7)
%   8
%   julia> nextpow(2, 9)
%   16
%   julia> nextpow(5, 20)
%   25
%   julia> nextpow(4, 16)
%   16
%   ```
%   """
assert(x>0 && a>1);
if x<=1
    outp = 1;
    return
end
n = ceil(log(x)/log(a));
p = a^(n-1);
if p>=x
    outp = p;
    return
else
    outp = a^n;
end
end