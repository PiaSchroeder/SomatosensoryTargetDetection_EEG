function [xp_perfam,r] = spm_dirichlet_exceedance_fam(alpha,Nsamp,partition)
% Compute exceedance probabilities for a Dirichlet distribution
% 
% Input:
% alpha     - Dirichlet parameters
% Nsamp     - number of samples used to compute xp [default = 1e6]
% partition - model partitioning, e.g. [1 1 1 2 3]
% 
% Output:
% xp_perfam - exceedance probability within model families
% r         - samples from posterior
%__________________________________________________________________________
%
% This function computes exceedance probabilities, i.e. for any given model
% k1, the probability that it is more likely than any other model k2.  
% More formally, for k1=1..Nk and for all k2~=k1, it returns p(x_k1>x_k2) 
% given that p(x)=dirichlet(alpha).
% This is done within model families, i.e. within each family, the
% exceedance probabilities sum up to unity.
% 
% Refs:
% Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ
% Bayesian Model Selection for Group Studies. NeuroImage (in press)
%__________________________________________________________________________
%
% Modified by PS: based on spm_dirichlet_exceedance.m, now computes XPs 
% within families and returns samples from the posterior
% distributions to compute family-level XPs in spm_BMS_family.m
%
% Original script: spm_dirichlet_exceedance.m
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% Will Penny & Klaas Enno Stephan
% $Id: spm_dirichlet_exceedance.m 3118 2009-05-12 17:37:32Z guillaume $


if nargin < 2
    Nsamp = 1e6;
end

Nk = length(alpha);                 % number of models
K  = length(unique(partition));     % number of families

% Size of families 
ind = cell(1,K);
fam_size = nan(1,K);
for i=1:K
    ind{i} = find(partition==i);
    fam_size(i) = length(ind{i});
end


% Sample from univariate gamma densities then normalise
% (see Dirichlet entry in Wikipedia or Ferguson (1973) Ann. Stat. 1,
% 209-230)
%----------------------------------------------------------------------
r = zeros(Nsamp,Nk);
for k = 1:Nk
    r(:,k) = spm_gamrnd(alpha(k),1,Nsamp,1);
end
sr = sum(r,2);
r = r./repmat(sr,1,Nk);

% Exceedance probabilities:
% For any given model k1, compute the probability that it is more
% likely than any other model k2~=k1
%----------------------------------------------------------------------
xp_perfam = nan(1,Nk);

% compute exceedance probabilities within families
for k = 1:K
    ri = r(:,ind{k});
    [~, j] = max(ri,[],2);
    xp = histc(j, 1:fam_size(k))';
    xp = xp / Nsamp;
    
    xp_perfam(1,ind{k}) = xp;
end
