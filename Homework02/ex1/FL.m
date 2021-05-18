function [ FL ] = FL(psi_i,psi_ii,u,i,N)

if i==1
    FL=(1/2)*(psi_i*(u(2)-u(1))-psi_ii*(u(1)-u(N)));
elseif i==N
    FL=(1/2)*(psi_i*(u(1)-u(N))-psi_ii*(u(N)-u(N-1)));
else
    FL=(1/2)*(psi_i*(u(i+1)-u(i))-psi_ii*(u(i)-u(i-1)));
end

end

