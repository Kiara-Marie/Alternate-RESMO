function out=knnp(ni,nf,II,minn,maxn,diffsn,T) % rate for transfer from n to n'
    
    % Rates calculated using PVS PRL 2008 
    % output units is um^3 ns^-1
    % use in conjunction with [ni,nf,II,minn,maxn,diffsn]=buildns(nl);

    kB = 1.3806504e-23;% #m2 kg s-2 K-1;J K-1; 8.62e-5 #eV K-1, % PVS PRL 2008
    Rydhc = 2.179872e-18;% #J: Ryd [cm-1] --> [J]

    eps_i=Rydhc./(power(ni,2.0)*kB*T);
    eps_f=Rydhc./(power(nf,2.0)*kB*T);
    max_eps_i=Rydhc./(power(minn,2)*kB*T);
    min_eps_i=Rydhc./(power(maxn,2)*kB*T);

    diffs=Rydhc.*diffsn./(kB*T); % scale dffsn properly

    out=II.*(kNOT(T).*power(eps_i,5/2).*power(eps_f,3/2)./power(max_eps_i,5/2))...
        .*exp(-(eps_i-min_eps_i)).*((22./(power(max_eps_i+0.9,7/3)))...
        +(4.5./(power(max_eps_i,2.5).*power(diffs+1-II,4/3)))); 

end    