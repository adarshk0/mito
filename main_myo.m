function ret = main_myo(tspan,y,XATPase)
% 	ics(1) = 9638.08;	% ATP
% 	ics(2) = 92.5336;	% ADP
% 	ics(3) = 6.61147;	% AMP
% 	ics(4) = 12000.1;	% PCr
% 	ics(5) = 12000;		% Cr 
% 	ics(6) = 2341.77;	% Pi 

	CONSTANTS = zeros(26,1);
	
	CONSTANTS(1) = 9670; %param
	CONSTANTS(2) = 9.0E2; %Kia
	CONSTANTS(3) = 3.49E4; %Kib
	CONSTANTS(4) = 1.55E4; %Kb
	CONSTANTS(5) = 2.224E2; %Kic
	CONSTANTS(6) = 3.49E4; %KIb
	CONSTANTS(7) = 4.73E3; %Kid
	CONSTANTS(8) = 1.67E3; %Kd
	CONSTANTS(9) = 13.00E3; %V1
	CONSTANTS(10) = 54.8E3; %V11
	CONSTANTS(11) = 63.8; %kfa
	CONSTANTS(12) = 1.68; %kba
	CONSTANTS(13) = 30; %t1
	CONSTANTS(14) = 60; %t2
	CONSTANTS(15) = 180; %t3
	CONSTANTS(16) = 100000; %ttt
	CONSTANTS(17) = 0.065800; %R
	CONSTANTS(18) = XATPase; %X_ATPase
	CONSTANTS(19) = 24.0; %KDT
	CONSTANTS(20) = 1.0E3; %Mgext
	CONSTANTS(21) = 347.0; %KDD
	CONSTANTS(22) = 380.0; %Mgx
	CONSTANTS(23) = (CONSTANTS(5).*CONSTANTS(8))./CONSTANTS(7); %Kc
	CONSTANTS(24) = 1; %G_CK
	CONSTANTS(25) = CONSTANTS(24); %G_AK
	CONSTANTS(26) = 1; %G_H
    
  	ics = y;
	
	% Solve model with ODE solver
 	[~,y] = ode15s(@(t,y) myo_odefun_mex(t,y,CONSTANTS), tspan, ics);
    ret = y(end,:);
end