function pde_out = pdesolve(t_in, STATES, L, extra)
	% Load stress
	yN = extra.yN;
	xN = extra.xN;
	uN = yN*xN;
    
    rshp = @(x) reshape(x,[yN,xN]);
    
    sol = StrangSplitting(t_in, STATES, L, extra);
    
    % Reshape output
	ATP		= rshp(sol((1:uN)+0*uN));
	ADP		= rshp(sol((1:uN)+1*uN));
	AMP		= rshp(sol((1:uN)+2*uN));
	PCr		= rshp(sol((1:uN)+3*uN));
	Cr		= rshp(sol((1:uN)+4*uN));
	Pi		= rshp(sol((1:uN)+5*uN));
	ATP_ig  = rshp(sol((1:uN)+6*uN));
	dPsi	= rshp(sol((1:uN)+7*uN));
	H_x		= rshp(sol((1:uN)+8*uN));
	NADH_x	= rshp(sol((1:uN)+9*uN));
	Pi_x	= rshp(sol((1:uN)+10*uN));
	QH2		= rshp(sol((1:uN)+11*uN));
	Cred	= rshp(sol((1:uN)+12*uN));
	O2		= rshp(sol((1:uN)+13*uN));
	K_x		= rshp(sol((1:uN)+14*uN));
	ATP_x	= rshp(sol((1:uN)+15*uN));
	ADP_x	= rshp(sol((1:uN)+16*uN));
    
	% Store outputs
	pde_out.ATP		= ATP;
	pde_out.ADP		= ADP;
	pde_out.AMP		= AMP;
	pde_out.PCr		= PCr;
	pde_out.Cr		= Cr;
	pde_out.Pi		= Pi;
	pde_out.ATP_ig	= ATP_ig;
	pde_out.dPsi	= dPsi;
	pde_out.H_x		= H_x;
	pde_out.NADH_x	= NADH_x;
	pde_out.Pi_x	= Pi_x;
	pde_out.QH2		= QH2;
	pde_out.Cred	= Cred;
	pde_out.O2		= O2;
	pde_out.K_x		= K_x;
	pde_out.ATP_x	= ATP_x;
	pde_out.ADP_x	= ADP_x;
    
end