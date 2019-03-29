function ktbr=kTBR(n,T)     % three-body-recombination
    
    % TBR rates output units in um^6 ns-1
    kB = 1.3806504e-23;% #m2 kg s-2 K-1;J K-1; 8.62e-5 #eV K-1,
    emass = 9.1093822e-31;% #kg
    h = 6.6260690e-34;% #m2 kg / s 4.13e-15 #eV s
    Rydhc = 2.179872e-18;% #J: Ryd [cm-1] --> [J]

    epsi=Rydhc./(power(n,2.0)*kB*T);
    lmbd=1e6*h./sqrt(2.0*pi*emass*kB*T);

    ktbr=kION(n,T).*power(n,2.0).*power(lmbd,3.0).*exp(epsi); % unit um^6 ns-1   *0.1 to test effect of TBR suppression
    ktbr(isnan(ktbr))=0;       % take care of computation error for values matlab sees too small and turns NaN
    ktbr(~isfinite(ktbr))=0;    % take care of computation error for values matlab sees too large and turns inf
end
