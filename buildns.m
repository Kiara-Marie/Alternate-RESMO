function [ni,nf,II,minn,maxn,diffsn]=buildns(nl)
    
    % Is needed to reevaluate the rate coeff.
    a=length(nl);
    II=ones(a,a); % will be matrix of ij=1 and ii=0
    ni=zeros(a,a);
    nf=zeros(a,a);
    minn=zeros(a,a);
    maxn=zeros(a,a);
    for i=1:a
        for j=1:a
            ni(i,j)=nl(i); % an array of initial states
            nf(i,j)=nl(j); % an array of final state
            minn(i,j)=min(ni(i,j),nf(i,j)); %find min of init and final state potential problem
            maxn(i,j)=max(ni(i,j),nf(i,j)); %find max of init and final state
            if i==j
            II(i,j)=0; % set to 0 ones with same init and final state
            end
        end
    end
    diffsn=abs(1./ni.^2-1./nf.^2); % difference in energy between the 2 states (no units)

end % help array