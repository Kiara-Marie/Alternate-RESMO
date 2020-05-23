function [T_PENNING,deac,D_DEAC_N_MIN,nden,NL,DEN0,NDEN,EDEN,DEAC]= penningIonization(ns,n_min,n0,den0,N,firstn,nl,Ry,kB,penningTimes)
% Sets initial conditions for electron and n-level distributions:
% no initial electrons -> calc. Penning seed electrons (look up Robicheaux)
penningByShell = zeros(N,1);
T_PENNING=zeros(1,N);
D_DEAC_N_MIN = zeros(n_min,N);
EDEN = zeros(N,1);
DEAC = zeros(N,1);
for ii=1:N %loop though all shells
    h=1+(ii-1)*ns;
    k=ii*ns;
    hh=1+(ii-1)*n_min;
    kk=ii*n_min;
   
    for jj = 1:penningTimes
        [PF,eden,rden]=penningfraction(n0,den0(ii));
        % Redistributes the Penning partners over lower n's:
        f=@(x)5.95*x.^5;                % This is the penning fraction distribution
        np=firstn:fix(n0/sqrt(2));      % Array of n states allowed after Penn ion
        ind=1:length(np);               % This is the distribution of penning fraction
        nden=nl*0;
        nden(ind)=eden*f(np/n0)/sum(f(np/n0));  % dividing by the sum normalizes the function
        nden(nl==n0)=rden;              % set n0 to rden
    end
    % Set initial temperature:    (Robicheaux 2005 JPhysB)
    T_PENNING(ii)=(-Ry*den0(ii)/n0^2 + Ry*rden/n0^2 + Ry*sum(nden(ind)./nl(ind).^2) )*1/(3/2*kB*eden); % by energy conservation
    deac=sum(nden(1:n_min));    % allow n<=n_min to decay
    D_DEAC_N_MIN(:,ii)=nden(1:n_min);
    nden(1:n_min)=zeros(n_min,1);
   
    NL(:,ii)=nl;       % save values for this shell in arrays
    DEN0(:,ii)=repmat(den0(ii),ns,1);
    NDEN(:,ii)=nden;
    EDEN(ii)=eden;
    DEAC(ii)=deac;
end
end