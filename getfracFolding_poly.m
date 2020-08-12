function [ff,pvals,Avals] = getfracFolding_poly(param)

% for a predefined set of parameters, calculate the steady-state
% concentrations
% and folding fraction
% by taking roots of a 4th degree polynomial
% ff(i) = fraction folding value for i-th valid solution
% pvals(i,:) = steady-state protein concentrations (Pg, Pgs, Pc, Pcs, P,
% Ps) for i-th valid solution
% Avals = all possible roots for alpha

% param is a structure which must have the following fields:
% k1, k2, k3, k4, k5, k6, k7, kd, kpgood, kps, c

% for folded proteins
m1 = -param.k1*param.k7/param.kd;
n1 = -param.k1-param.k2-param.k6;
n2 = param.k6 - param.k1*param.k4/param.kd;
b1 = -param.k7*param.kpgood/param.kd;
b2 = -param.kpgood*(1+param.k4/param.kd);
% determinant polynomial
detcoeff = [m1*param.k3, n2*param.k3 + n1*param.k3+m1*param.k5, n1*param.k5];
% solving for Pc
%factor = (1+param.k1/param.k2);
factor=1;
pccoeff = [factor*b1*param.k3, factor*b1*param.k5+b2*param.k3,0];
% solving for Pg
pgcoeff = [n2*b1-m1*b2, -n1*b2];



%checkpc = polyval(pccoeff,alpha)/polyval(detcoeff,alpha)

% for unfoldable proteins
dE = param.de;

% paramS.k2 = param.k2*exp(dE);
% paramS.k3 = param.k3*exp(-dE);
% paramS.k4 = param.k4;
% paramS.k5 = param.k5;
% paramS.k6 = param.k6;
% paramS.k7 = param.k7;

paramS.k2 = param.k2*exp(param.de);
paramS.k3 = param.k3;
paramS.k4 = param.k4*exp(-param.de);
paramS.k5 = param.k5;
paramS.k6 = param.k6;
paramS.k7 = param.k7;

n1s = -(paramS.k2+paramS.k6);
b1s = -paramS.k7*param.kps/param.kd;
b2s = -param.kps*(1+paramS.k4/param.kd);
detScoeff = [paramS.k3*paramS.k6 + n1s*paramS.k3, n1s*paramS.k5];
pcScoeff = [b1s*paramS.k3, b1s*paramS.k5 + b2s*paramS.k3,0];
pgScoeff = [b1s*paramS.k6, -n1s*b2s];
%checkpcS = polyval(pcScoeff,alpha)/polyval(detScoeff,alpha);

%% automatically set up quartic polynomial
% multiply polynomials by convolution
detdet = conv(detcoeff,detScoeff); % product of determinant polynomials
detdetA = conv(detdet,[1+paramS.k7*param.Pb/paramS.k2 0]);

%this factor accounts for a folded protein needing a process at
%rate k2 to detach from the chaperone
factor = (1+param.k1/(param.k2+param.k6));
qd = factor*[0 conv(pccoeff,detScoeff)];
qdS = conv(pcScoeff,detcoeff);

% final polynomial for alpha
polyA = detdetA +qd + qdS - param.c*[0,detdet];

%polyA

% possible values for alpha
Avals = roots(polyA);
%Avals
% allowable roots for alpha
ind = find(Avals > 0 & Avals <=param.c);

ff = []; pvals = [];
ct = 0;
for ac = ind
    % get steady state protein concentrations
    alpha = Avals(ac);
    pc = polyval(pccoeff,alpha)/polyval(detcoeff,alpha);
    pcs = polyval(pcScoeff,alpha)/polyval(detScoeff,alpha);
    pg = polyval(pgcoeff,alpha)/polyval(detcoeff,alpha);
    pgs = polyval(pgScoeff,alpha)/polyval(detScoeff,alpha);
    p = (param.kpgood - param.k1*pc)/param.kd;
    ps = param.kps/param.kd;
    plist = [pg pgs pc pcs p ps];
    
    if (all(plist>0) & pc < param.c & pcs < param.c)
        % valid set of parameters
        ct = ct+1;
        pvals(ct,:) = plist;
        ff(ct) = param.k1*pc/param.kpgood;
    end
end
