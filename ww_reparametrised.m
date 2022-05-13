function ww = ww_reparametrised(zz, par, inverse)
if nargin<3
    inverse=0;
end
if ~inverse
    switch par.repar
        case 'pos'
            ww = exp(zz);
        case 'possq'
            ww = zz.^2;
        case 'posb'
            ww = par.wmax*exp(-exp(zz));
        otherwise
            ww=zz;
    end
else
    switch par.repar
        case 'pos'
            ww = log(zz);
        case 'possq'
            ww = sqrt(zz);
        case 'posb'
            ww = log(-log(zz/par.wmax));
        otherwise
            ww=zz;
    end
end
end