function [T_PENNING,deac,D_DEAC_N_MIN,nden,NL,DEN0,NDEN,EDEN,DEAC]= penningIonization(ns,n_min,n0,den0,N,firstn,nl,Ry,kB)
% Sets initial conditions for electron and n-level distributions:
% no initial electrons -> calc. Penning seed electrons (look up Robicheaux)
T_PENNING=zeros(1,N);
D_DEAC_N_MIN = zeros(n_min,N);
EDEN = zeros(N,1);
DEAC = zeros(N,1);
NDEN = zeros(length(nl),N);
DEAC = zeros(N,1);
NL = zeros(length(nl), N);
for ii=1:N %loop though all shells
    h=1+(ii-1)*ns;
    k=ii*ns;
    hh=1+(ii-1)*n_min;
    kk=ii*n_min;
    nden=nl*0;
    f=@(x)5.95*x.^5;                % This is the penning fraction distribution
    
    [PF,primaryEden,primaryRden]=penningfraction(n0,den0(ii));
    % Redistributes the Penning partners over lower n's:  
    np=firstn:fix(n0/sqrt(2));      % Array of n states allowed after Penn ion
    ind=1:length(np);               % This is the distribution of penning fraction
    % dividing by the sum normalizes the function,
    % multiplying by eden gives it the correct magnitude
    nden(ind) = nden(ind) + (primaryEden*f(np/n0)/sum(f(np/n0)))';
    nden(nl==n0)= primaryRden;              % set n0 to rden
    nden(1:n_min)=zeros(n_min,1);
%     figure(1);
%     ax = gca;
%     histogram(ax,'BinEdges',nl(1:n0+1)','BinCounts',nden(1:n0));
%     hold on;
    spMean = round(sum(nden(ind).*np')/sum(nden(ind)));
    
    
    penningPartnerProportion = primaryEden/(primaryEden+primaryRden);
    remainingProportion = primaryRden/(primaryEden+primaryRden);
    n_avg = round((spMean * penningPartnerProportion) + (n0 * remainingProportion));
    combinedDen = den0(ii) - primaryEden;
    [PF, secondaryEden, secondaryRden] = penningfraction(n_avg, combinedDen);
    
    % multiply by secondaryEden, scaled by which proportion were n0 gives it the correct magnitude
    nden(ind)= nden(ind) + (secondaryEden*remainingProportion*f(np/n0)/sum(f(np/n0)))';
    nden(ind)= nden(ind) - (secondaryEden*penningPartnerProportion*f(np/n0)/sum(f(np/n0)))';
    nden(nl==n0) = nden(nl==n0) - secondaryEden*remainingProportion;
    nden(1:n_min)=zeros(n_min,1);
    
    secondaryPenningN=firstn:fix(spMean/sqrt(2));      % Array of n states allowed after Penn ion
    spInd = 1:length(secondaryPenningN);
    secondaryPenningDist = f(secondaryPenningN/spMean)/sum(f(secondaryPenningN/spMean));
    nden(spInd) = nden(spInd) + (secondaryEden*penningPartnerProportion*secondaryPenningDist)';
    
%     histogram(ax,'BinEdges',nl(1:n0+1)','BinCounts',nden(1:n0));    
%     eden = primaryEden + secondaryEden;
    
    % Set initial temperature:    (Robicheaux 2005 JPhysB)
    T_PENNING(ii)=(-Ry*den0(ii)/n0^2 + Ry*nden(n0)/n0^2 + Ry*sum(nden(nl < n0)./nl(nl < n0).^2) )*1/(3/2*kB*eden); % by energy conservation
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