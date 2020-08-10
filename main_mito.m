function sol = main_mito(tspan,y, M_basal, mass, isOnBoundary)
% Initialise constants and state variables
%	ics(1) = 9610.35;	% ATP
%	ics(2) = 113.792;	% ADP
%	ics(3) = 13.0731;	% AMP
%	ics(4) = 11998.6;	% PCr
%	ics(5) = 12001.4;	% Cr
%	ics(6) = 2306.24;	% Pi
%	ics(7) = 200007.8;	% ATP_ig
%	ics(8) = 174.993;	% dPsi
%	ics(9) = 0.0725301; % H_x
%	ics(10) = 1552.31;	% NADH_x
%	ics(11) = 2058.61;	% Pi_x
%	ics(12) = 597.711;	% QH2
%	ics(13) = 439.884;	% Cred
%	ics(14) = 47.25;	% O2
%	ics(15) = 137085;	% K_x
%	ics(16) = 8964.91;	% ATP_x
%	ics(17) = 1035.09;	% ADP_x

CONSTANTS = zeros(78,1);

ics = y;

CONSTANTS(1) = 6.6E3; % V1
CONSTANTS(2) = 30.41E3; % V11
CONSTANTS(3) = 7.5E2; % Kia
CONSTANTS(4) = 2.05E2; % Kic
CONSTANTS(5) = 5.0E2; % Kd
CONSTANTS(6) = 5.2E3; % Kb
CONSTANTS(7) = 2.6E4; % Kib
CONSTANTS(8) = 1.6E3; % Kid
CONSTANTS(9) = 2.6E4; % KIb
CONSTANTS(10) = 17; % KDTi
CONSTANTS(11) = 380.0; % Mgi
CONSTANTS(12) = 282.0; % KDDi
CONSTANTS(13) = 17.0; % KDTi
CONSTANTS(14) = 0.0083; % R
CONSTANTS(15) = 298.0; % T
CONSTANTS(16) = 0.0965; % F
CONSTANTS(17) = 0.861; % u
CONSTANTS(18) = 10E6; % eta
CONSTANTS(19) = 2.5; % na
CONSTANTS(20) = 20; % param
CONSTANTS(21) = 2.4734; % RT
CONSTANTS(22) = 0.096484; % F
CONSTANTS(23) = 3; % n_A
CONSTANTS(24) = -69.37; % dG_C1o
CONSTANTS(25) = -32.53; % dG_C3o
CONSTANTS(26) = -122.94; % dG_C4o
CONSTANTS(27) = 36.03; % dG_F1o
CONSTANTS(28) = 7.1; % pH_e
CONSTANTS(29) = 0.15e6; % K_e
CONSTANTS(30) = 0; % ATP_e
CONSTANTS(31) = 0; % AMP_e
CONSTANTS(32) = 2.4e1; % K_DT
CONSTANTS(33) = 3.47e2; % K_DD
CONSTANTS(34) = 0.4331; % K_AK
CONSTANTS(35) = 0.72376; % W_m
CONSTANTS(36) = 5.99; % gamma
CONSTANTS(37) = 0.0027e6; % Ctot
CONSTANTS(38) = 0.00135e6; % Qtot
CONSTANTS(39) = 0.00297e6; % NADtot
CONSTANTS(40) = 1.3413e2; % k_Pi1
CONSTANTS(41) = 6.7668e2; % k_Pi2
CONSTANTS(42) = 1.9172e2; % k_Pi3
CONSTANTS(43) = 0.02531e6; % k_Pi4
CONSTANTS(44) = 4.5082e2; % k_PiH
CONSTANTS(45) = 4.5807; % r
CONSTANTS(46) = 0.09183; % x_DH
CONSTANTS(47) = 0.36923; % x_C1
CONSTANTS(48) = 0.091737; % x_C3
CONSTANTS(49) = 3.2562e-5; % x_C4
CONSTANTS(50) = 150.93e-6; % x_F1
CONSTANTS(51) = 12.0e3; % x_ANT
CONSTANTS(52) = 339430; % x_Pi1
CONSTANTS(53) = 2.9802e1; % x_KH
CONSTANTS(54) = 250; % x_Hle
CONSTANTS(55) = 0; % x_K
CONSTANTS(56) = 3.5e0; % k_mADP_i
CONSTANTS(57) = 0; % x_AK
CONSTANTS(58) = 85; % p_A
CONSTANTS(59) = 1.2e2; % k_O2
CONSTANTS(60) = 100e-6; % x_buff
CONSTANTS(61) = 327; % x_Pi2
CONSTANTS(62) = 1e-6; % mincond
CONSTANTS(63) = 63.8; % kfa
CONSTANTS(64) = 1.68; % kba
CONSTANTS(65) = 6.756756756756757; % C_im
CONSTANTS(66) = 17.0; % KDTx
CONSTANTS(67) = 380.0; % Mgx
CONSTANTS(68) = 282.0; % KDDx
CONSTANTS(69) = 1.0E3; % Mgext
CONSTANTS(70) = 0.01; % Rexch
CONSTANTS(71) = 10.^(6-CONSTANTS(28)); % H_e
CONSTANTS(72) = 10^(-0.48); % k_dHatp
CONSTANTS(73) = 10^(-0.29); % k_dHadp
CONSTANTS(74) = CONSTANTS(29); % K_i
CONSTANTS(75) = 0.9.*CONSTANTS(35); % W_x
CONSTANTS(76) = 0.1.*CONSTANTS(35); % W_i
CONSTANTS(77) = CONSTANTS(71); % H_i
CONSTANTS(78) = 10^(-0.75); % k_dHPi

% Apply feedback	
CONSTANTS(47) = (mass/M_basal)*CONSTANTS(47); % Modify x_C1
CONSTANTS(48) = (mass/M_basal)*CONSTANTS(48); % Modify x_C3
CONSTANTS(49) = (mass/M_basal)*CONSTANTS(49); % Modify x_C4
CONSTANTS(50) = (mass/M_basal)*CONSTANTS(50); % Modify x_F1
CONSTANTS(54) = (2*M_basal/(M_basal+mass))*CONSTANTS(54); % Modify X_Hle (X_leak)

% Set numerical accuracy options for ODE solver
% options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1);

% Solve model with ODE solver
[~,y] = ode15s(@(t,y) mito_odefun(t, y, CONSTANTS, isOnBoundary), tspan, ics);
sol = y(end,:);

end
