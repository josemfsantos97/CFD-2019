function [ psi ] = psi_minmod( r )

if isnan(r)
    psi=0;
elseif r<0
    psi=0;
elseif r>=0 && r<=1
    psi=r;
elseif r>1
    psi=1;
end

end

