function [vp,vm]=fluxsplitting_scalar_general(u,flux, dflux,strategy)
% split the flux into right-going and left-going
% OUTPUT:
%   * vp: positive dflux, v^{+}, which corresponds to f_{i+1/2}^{-}
%   * vn: negative dflux  v^{-}, which corresponds to f_{i+1/2}^{+}

a=max(abs(dflux(u)));

switch strategy
    case 'Upwind'
        vp= max(u,0); %dflux^{+}
        vm=-min(u,0); %dflux^{-}
    case 'LF'
        vp =  0.5*(flux(u)+a*u); %dflux^{+}
        vm = -0.5*(flux(u)-a*u); %dflux^{-}
    otherwise
        error('only case not available')
end