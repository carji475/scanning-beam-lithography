function [cost,grad] = exposureFunction_lbfgs_fft(varargin)
%{
    Calculation of cost and gradient using the FFT
%}
persistent a tr Nx Ny scale lambda gamma par
if nargin>3
    a      =varargin{1};
    tr     =varargin{2};
    Nx     =varargin{3};
    Ny     =varargin{4};
    scale  =varargin{5};
    lambda =varargin{6};
    gamma  =varargin{7};
    par    =varargin{8};
else
    zz = varargin{1};
    Z  = varargin{2};
    indices = varargin{3};
    ww = ww_reparametrised(zz, par);
    
    XX = fft_exposure(squareUp(ww,indices,Nx,Ny));
    Zh = 1./(1+exp(-a*(XX-tr)));
    pe = Z(:)-Zh(:);
    
    % cost function
    cost = scale*(pe'*pe + lambda*sum(ww(:)) + gamma*(XX(:)'*XX(:)));
    
    % gradient
    if nargout>1
        grad = fft_exposure(-2*a*(Z-Zh).*Zh.*(1-Zh) + 2*gamma*XX)*scale;
        grad = grad(indices)+lambda*ones(sum(indices),1);
        switch par.repar
            case 'pos'
                grad=grad.*ww;
            case 'possq'
                grad=grad*2.*zz;
            case 'posb'
                grad=grad.*(-ww.*exp(zz));
            otherwise
        end
    end
end
end
