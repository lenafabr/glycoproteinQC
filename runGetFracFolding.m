clear all;
close all;

param.k1 = 1; %this is k_f
param.k2 = 1; %this is k_r
param.k3 = 1; %this is k_c
param.k4 = 1; %this is k_g
param.k5 = 1; %this is k_-g
param.k6 = 1; %this is k_-c
param.k7 = 1; %this is k_-r
param.kd = 1; %this is k_d
param.kpgood = 1; %this is k_p
param.kps = 1; % this is k_p*
param.c = 1; %this is C_tot
param.Pb = 0; %this is P_b
%this is a parameter that is not used in the manuscript, it should be kept
% at zero
param.de = 0; 

[ff,pvals,Avals] = getfracFolding_poly(param);
%ff is the folding fraction

%pvals is a vector of protein concentrations
%[Pg Pg* Pc Pc* P P*]

%Avals is a vector of the available chaperone concentration C_A
%corresponding to the roots of the solved quartic equation. only the valid
%(C_A > 0 and C_A < C_tot) solution is used