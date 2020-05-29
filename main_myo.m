function ret = main_myo(tspan,y,XATPase)
	CONSTANTS = zeros(26,1);
    
	CONSTANTS(1) = 9670; % param	
	CONSTANTS(2) = 9.0E2; % Kia		
	CONSTANTS(3) = 3.49E4; % Kib		
	CONSTANTS(4) = 1.55E4; % Kb		
	CONSTANTS(5) = 2.224E2; % Kic		
	CONSTANTS(6) = 3.49E4; % KIb		
	CONSTANTS(7) = 4.73E3; % Kid		
	CONSTANTS(8) = 1.67E3; % Kd		
	CONSTANTS(9) = 13.00E3; % V1
	CONSTANTS(10) = 54.8E3; % V11		
	CONSTANTS(11) = 63.8; % kfa		
	CONSTANTS(12) = 1.68; % kba		
	CONSTANTS(13) = 0.03; % t1		
	CONSTANTS(14) = 0.06; % t2		
	CONSTANTS(15) = 0.18; % t3		
	CONSTANTS(16) = 100; % ttt		
	CONSTANTS(17) = 0.065800; % R		
	CONSTANTS(18) = XATPase*1e6; %0.05*1e6; % X_ATPase	
	CONSTANTS(19) = 24.0; % KDT		
	CONSTANTS(20) = 1.0E3; % Mgext	
	CONSTANTS(21) = 347.0; % KDD		
	CONSTANTS(22) = 380.0; % Mgx
	CONSTANTS(23) = 1; % G_CK
	CONSTANTS(24) = ( CONSTANTS(5).*CONSTANTS(8))./CONSTANTS(7); % Kc
	CONSTANTS(25) = CONSTANTS(23); % G_AK		
	CONSTANTS(26) = 1; % G_H	

  	ics = y;
	
	% Solve model with ODE solver
 	[~,y] = ode15s(@(t,y) myo_odefun_mex(t,y,CONSTANTS), tspan, ics);
    ret = y(end,:);
end