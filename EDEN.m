function [ ed, x ] = EDEN(eden, vol, shell, time)
%%%%%%%%%%%%%% use this function as follows
% [y, x]=EDEN(eden, vol, rx, 1);    plot density profile over all shells at t=1;
% [y]=EDEN(eden, vol);              plot total electron number over time;
% [y]=EDEN(eden, vol, 1);           plot total electron number in shell 1 over time;


%   eden(time,shell)
res= vol.*eden;
if nargin < 3
    ed= sum(res,2);
elseif nargin < 4
    ed=res(:,shell);
else
    edr=eden(time,:);
    ed=[flip(edr) edr];
    xr=shell(time,:);
    x=[-flip(xr) xr];
end

end
