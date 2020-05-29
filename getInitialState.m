function [t_index, fuse_event, split_event, biogen_event, stressor_event, mitochondria, damage, atp, adp, amp, pcr, cr, pi, dpsi, h_x, nadh_x, pi_x, qh2, cred, o2, k_x, atp_x, adp_x] = getInitialState(input0,input)

sz = [21, 79];
xN = sz(2);
yN = sz(1);
uN = xN*yN;

dx = 1;
dy = 0.75;
		
dt = 0.01;
tf = 120;
tmesh = 0:dt:tf;
tN = length(tmesh);

% Declare and initialise state variables
Z = zeros(yN,xN,tN);
[mitochondria, atp, adp, amp, pcr, cr, pi, dpsi, h_x, nadh_x, pi_x, qh2, cred, o2, k_x, atp_x, adp_x] = deal(Z);
damage = NaN(yN,xN,tN);
% Keep track of ABM events
fuse_event     = zeros(1,tN);
split_event    = zeros(1,tN);
biogen_event   = zeros(1,tN);
stressor_event = zeros(1,tN);

if nargin>1
	t_index = size(input.M,3)+1;
	for t=1:size(input.M,3)
		fuse_event(t)		= input.fuse(t);
		split_event(t)		= input.split(t);
		biogen_event(t)		= input.biogen(t);
		stressor_event(t)	= input.stress(t);
		mitochondria(:,:,t)	= input.M(:,:,t);
		damage(:,:,t)		= input.d(:,:,t);
		atp(:,:,t)	 		= input.ATP(:,:,t);
		adp(:,:,t)	 		= input.ADP(:,:,t);
		amp(:,:,t)	 		= input.AMP(:,:,t);
		pcr(:,:,t)	 		= input.PCr(:,:,t);
		cr(:,:,t)	 		= input.Cr(:,:,t);
		pi(:,:,t)	 		= input.Pi(:,:,t);
		dpsi(:,:,t)	 		= input.dPsi(:,:,t);
		h_x(:,:,t)	 		= input.H_x(:,:,t);
		nadh_x(:,:,t)		= input.NADH_x(:,:,t);
		pi_x(:,:,t)	 		= input.Pi_x(:,:,t);
		qh2(:,:,t)	 		= input.QH2(:,:,t);
		cred(:,:,t)	 		= input.Cred(:,:,t);
		o2(:,:,t)	 		= input.O2(:,:,t);
		k_x(:,:,t) 	 		= input.K_x(:,:,t);
		atp_x(:,:,t) 		= input.ATP_x(:,:,t);
		adp_x(:,:,t)		= input.ADP_x(:,:,t);
	end
else
	t_index = 2;
	M0 = dlmread('ic.txt');
	% Declare and initialise state variables
	Z = zeros(yN,xN,tN);
	im = logical(M0);

	mitochondria = Z;
	mitochondria(:,:,1) = M0;

	[atp, adp, amp, pcr, cr, pi, dpsi, h_x, nadh_x, pi_x, qh2, cred, o2, k_x, atp_x, adp_x] = deal(Z);

	damage_init = input0;

	MYO_INIT_STATES(1) = 9638.08; % ATP
	MYO_INIT_STATES(2) = 92.5336; % ADP
	MYO_INIT_STATES(3) = 6.61147; % AMP
	MYO_INIT_STATES(4) = 12000.1; % PCr
	MYO_INIT_STATES(5) = 12000; % Cr
	MYO_INIT_STATES(6) = 2341.77; % Pi
	MYO_INIT_STATES(7) = 47.25; % O2

	MIT_INIT_STATES(1) = 9607.15; % ATPi
	MIT_INIT_STATES(2) = 112.217; % ADPi
	MIT_INIT_STATES(3) = 17.8491; % AMPi
	MIT_INIT_STATES(4) = 12000; % PCri
	MIT_INIT_STATES(5) = 12000; % Cri
	MIT_INIT_STATES(6) = 2328.17; % Pii
	MIT_INIT_STATES(7) = 167.954; % dPsi
	MIT_INIT_STATES(8) = 0.0724323; % H_x
	MIT_INIT_STATES(9) 	= 0.00154986e6; % NADH_x
	MIT_INIT_STATES(10) = 0.0020e6; % Pi_x
	MIT_INIT_STATES(11) = 0.000599417e6; % QH2
	MIT_INIT_STATES(12) = 0.000438911e6; % Cred
	MIT_INIT_STATES(13) = 47.25; % O2
	MIT_INIT_STATES(14) = 0.137087e6; % K_x
	MIT_INIT_STATES(15) = 0.009e6; % ATP_x
	MIT_INIT_STATES(16) = 0.001e6; % ADP_x


	damage = NaN(yN,xN,tN);
	d0 = damage(:,:,1);
	d0(logical(M0)) = 0;
	d0(2:6,:) = d0(2:6,:) + damage_init;
	damage(:,:,1) = d0;

	atp = init_consts(atp, MIT_INIT_STATES(1), MYO_INIT_STATES(1), im);
	adp = init_consts(adp, MIT_INIT_STATES(2), MYO_INIT_STATES(2), im);
	amp = init_consts(amp, MIT_INIT_STATES(3), MYO_INIT_STATES(3), im);
	pcr = init_consts(pcr, MIT_INIT_STATES(4), MYO_INIT_STATES(4), im);
	cr  = init_consts(cr,  MIT_INIT_STATES(5), MYO_INIT_STATES(5), im);
	pi  = init_consts(pi,  MIT_INIT_STATES(6), MYO_INIT_STATES(6), im);
	o2	= init_consts(o2,  MIT_INIT_STATES(13),MYO_INIT_STATES(7), im);

	dpsi 	= init_consts(dpsi, 	MIT_INIT_STATES(7),  0, im);
	h_x	 	= init_consts(h_x,  	MIT_INIT_STATES(8),  0, im);
	nadh_x  = init_consts(nadh_x, 	MIT_INIT_STATES(9),  0, im);
	pi_x 	= init_consts(pi_x, 	MIT_INIT_STATES(10), 0, im);
	qh2	 	= init_consts(qh2,  	MIT_INIT_STATES(11), 0, im);
	cred 	= init_consts(cred, 	MIT_INIT_STATES(12), 0, im);
	k_x	 	= init_consts(k_x, 		MIT_INIT_STATES(14), 0, im);
	atp_x 	= init_consts(atp_x,	MIT_INIT_STATES(15), 0, im);
	adp_x 	= init_consts(adp_x,	MIT_INIT_STATES(16), 0, im);

	tmp=[];
	% Keep track of ABM events
	fuse_event     = zeros(1,tN);
	split_event    = zeros(1,tN);
	biogen_event   = zeros(1,tN);
	stressor_event = zeros(1,tN);
end
