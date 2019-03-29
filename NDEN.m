function [ nd, x ] = NDEN(nden, vol, shell, time)
%%%%%%%%%%%%%% use this function as follows
% [y, x]=NDEN(nden, vol, rx, 1);    plot density profile over all shells at t=1;
% [y]=NDEN(nden, vol);              plot total electron number over time;
% [y]=NDEN(nden, vol, 1);           plot total electron number in shell 1 over time;

timesteps=size(vol,1);
N=size(vol,2); %get number of shells from vol array size
ns=size(nden,2)/N; %number of n-levels (usually n=1,...,100)
VOL=zeros(size(nden));
for ii=1:N
    VOL(:,1+(ii-1)*ns:ii*ns)=bsxfun(@times,ones(timesteps,ns),vol(:,ii));%vol(:,ii).*ones(timesteps,ns);  %build VOL matrix
end

%   eden(time,shell)
res= VOL.*nden;
if nargin < 3
    nd= sum(res,2);
elseif nargin < 4
    nd=sum(res(:,ns*(shell-1)+1:ns*shell),2);
else
    ndr=zeros(1,N);
    for ii=1:N
        ndr(ii)=sum(nden(time,1+(ii-1)*ns:ii*ns),2);
    end
    nd=[flip(ndr) ndr];
    xr=shell(time,:); %now shell contains the positions
    x=[-flip(xr) xr];
end

end
