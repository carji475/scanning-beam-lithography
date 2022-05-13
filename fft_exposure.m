function conv = fft_exposure(varargin)
persistent M1 M2 L1 L2 N1 N2 P1 P2 kzfft
if nargin>1
    M1  = varargin{1};
    M2  = varargin{2};
    support_x  = varargin{3};
    support_y  = varargin{4};
    xres = varargin{5};
    yres = varargin{6};
    Hsqn = varargin{7};
    L1 = min(M1-1, support_x );
    L2 = min(M2-1, support_y );
    N1 = 2*L1+1;
    N2 = 2*L2+1;
    P1 = nextprod([2 3 5 7], M1+N1);
    P2 = nextprod([2 3 5 7], M2+N2);
    kz = exp( -0.5*( reshape(sum( (Hsqn*[repmat( (-L1:L1)*xres, 1, N2 ); ...
        repelem( (-L2:L2)*yres, 1, N1 )]).^2 ), N1, N2) ) ); % Gaussian power kernel
    kzfft = fft2(kz,P1,P2);
else
    cc = zeros(P1,P2);
    cc(L1+1:L1+M1,L2+1:L2+M2) = varargin{1};
    
    conv = ifft2( fft2(cc).*kzfft );
    conv = conv(N1:N1+M1-1,N2:N2+M2-1);
end
end