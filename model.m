function model_output = model
% pars = [fuse, split, XATPase, O2, damage]
index_input = 1; % filename 
par = [1,1,0.05,47.25,0]; % model parameters
mach_tol = 1e-9;

% Generate spatial mesh:
% 	Our initial geometry is from experimental images
%	Our discretisation is based on the average size of a mito
dx = 1;
dy = 0.75;

xN = 79;
yN = 21;
uN = xN*yN;

dt = 0.01;
tf = dt;
% Generate time mesh:
% 	Simulate our model over the period of [0, tf]
% 	The average duration of murine cardiac cycle is 180 milliseconds
tmesh = 0:dt:tf;
tN = length(tmesh);

% Construct diffusion operators
% Dirichlet BC
Ex=sparse(2:xN,1:xN-1,1,xN,xN);
Ax=Ex+Ex'-2*speye(xN);
Ax(1,1)=-1; Ax(xN,xN)=-1;  %Dirichlet B.Cs
Ey=sparse(2:yN,1:yN-1,1,yN,yN);
Ay=Ey+Ey'-2*speye(yN);
Ay(1,1)=-1; Ay(yN,yN)=-1; %Dirichlet B.Cs
L_diri=kron(Ax/dx^2,speye(yN))+kron(speye(xN),Ay/dy^2);
bc=sparse(yN,xN);
bc(1,:) = 1; bc(end,:) = 1; 
bc(:,1)=1; bc(:,end) =1;
ix = bc==1;
L_diri(ix,:) = 0;

% No flux
Ex=sparse(2:xN,1:xN-1,1,xN,xN);
Ax=Ex+Ex'-2*speye(xN);
Ax(1,1)=-1; Ax(xN,xN)=-1;  %Neumann B.Cs
Ey=sparse(2:yN,1:yN-1,1,yN,yN);
Ay=Ey+Ey'-2*speye(yN);
Ay(1,1)=-1; Ay(yN,yN)=-1; %Neumann B.Cs
L_nf=kron(Ax/dx^2,speye(yN))+kron(speye(xN),Ay/dy^2);

Zh = sparse(size(L_nf,1), size(L_nf,2));

% Diffusion coefficients
diff_ANP = 30;
diff_cr  = 260;
diff_pi  = 327;
diff_o   = 300;

% Load ABM parameters
M_single = 0.825;
M_max = 43;
M_basal = 22;

lam_fuse = 0.0167*par(1);
lam_split = 2e-3*par(2);
lam_bio = 5.7670e-04;
k_damage = 1e-3;

damage_init = par(5);
Td = 10;
p_damage = 0.01;
p_death  = 0.6;
dam_c = 250;

IE0 = 0.05;

% Workload (VO2)
XATPase = par(3);

% Oxygen boundary
o2_val = par(4);

% Get inital states
[t_index, fuse_event, split_event, biogen_event, stressor_event, mitochondria, damage,...
 atp, adp, amp, pcr, cr, pi, atp_ig, dpsi, h_x, nadh_x, pi_x, qh2, cred, o2, k_x, atp_x, ...
 adp_x] = getInitialState([damage_init,tf]);
 
% Imposing bdry conditions
%o2 = o2/47.25;
%o2 = o2_val*o2;

o2(1,:) = o2_val;
o2(end,:) = o2_val;
o2(:,1) = o2_val;
o2(:,end) = o2_val;

% Evolve our ABM
for t=t_index:tN
    A 		= sparse(mitochondria(:,:,t-1));
    dam		= damage(:,:,t-1);
    ATP		= atp(:,:,t-1);
    ADP		= adp(:,:,t-1);
    AMP		= amp(:,:,t-1);
    PCr		= pcr(:,:,t-1);
    Cr		= cr(:,:,t-1);
    Pi		= pi(:,:,t-1);
    ATP_ig	= atp_ig(:,:,t-1);
    dPsi	= dpsi(:,:,t-1);
    H_x     = h_x(:,:,t-1);
    NADH_x	= nadh_x(:,:,t-1);
    Pi_x	= pi_x(:,:,t-1);
    QH2		= qh2(:,:,t-1);
    Cred	= cred(:,:,t-1);
    O2		= o2(:,:,t-1);
    K_x		= k_x(:,:,t-1);
    ATP_x	= atp_x(:,:,t-1);
    ADP_x	= adp_x(:,:,t-1);
    
       
    % Fusion events:
    % 	Each event requires the energetic stress to be recalculated
    %	Those agents who will fuse are selected based on a prob.
    %	These agents are then iterated over
    ids = unique(A(A>0)); % generate ID's for agents
    Es = zeros(1, numel(ids));
    pFuse =  zeros(1, numel(ids));
    for n=1:numel(ids)
        Es(n) = get_stress(ids(n));
        pFuse(n) = (1-exp(-(lam_fuse+Es(n)./(IE0+Es(n)))*dt));
    end
    stressor_event(t) = sum(Es); % log event

    will_fuse = rand(size(pFuse))<=pFuse;
    idf = ids(will_fuse); % ids of those that fuse
    fuse_event(t) =sum(will_fuse); % log the number of fusion events 		
    ind = 1;
    % Iterate over chosen mitochondria
    while numel(idf)>0 && ind<=numel(idf)
        n = idf(ind);
        % remove the index
        idf(idf==n) =[];

        % get adjacent mitochondria
        adj_mito = reshape(get_adj(n),1,[]);
        to_fuse = adj_mito;
        % if there's no adjacent mitochondria then skip to the next chosen mito
        if ~numel(adj_mito)
            continue
        end

        % check the masses are less than M_max otherwise remove from list
        for cnter=1:numel(adj_mito)
            if get_mass(n)+get_mass(adj_mito(cnter))>=M_max
                to_fuse(to_fuse==adj_mito(cnter))=[];
            end
        end

        % Updated damage is the average damage but also cut off where appropriate
        new_dam = dam(ismember(A,to_fuse));
        if mean(new_dam(:))<1
            dam(ismember(A,to_fuse)) = 0;
            dam(ismember(A,n)) =0;
        else
            dam(ismember(A,to_fuse)) = mean(new_dam(:));
            dam(ismember(A,n))= mean(new_dam(:));
        end
        A(ismember(A,to_fuse)) = n;

        % remove any mitos that need to fuse that happen to be ones we've already fused
        idf(ismember(idf,to_fuse)) = [];

        %update counter
        ind=ind+1;

    end
	abc = 1;
    % Fission events
    ids = unique(A(A>0)); % get different ids for mitos
    Es = zeros(1, numel(ids));
    pSplit =  zeros(1, numel(ids));
    for n=1:numel(ids)
        Es(n) = get_stress(ids(n));
        kdam = get_damage(ids(n));
        Mass = get_mass(ids(n));
        pSplit(n) = max(1-exp(-lam_split*Mass*dt),kdam./(dam_c+kdam));
    end

    will_split = rand(size(pSplit))<=pSplit;
    idf = ids(will_split); % ids of those that split
    split_event(t) =sum(will_split); % log splitting events
    ind = 1;
    while numel(idf)>0 && ind<=numel(idf)
        n = idf(ind);

        % remove the index from the list
        idf(idf==n) =[];

        % if the mass is less than twice the min then skip
        if get_mass(n)<= 2*M_single
            continue
        end

        % the new id of the mitochondria
        new_id = get_new_id(A);

        % calculate the new mass of the mito
        M1  = M_single + (get_mass(n)-2*M_single)*rand;

        % actually split the mito
        idx = find(A==n);
        A(idx(1:round(M1/M_single))) = new_id;

        % update damage
        es = get_stress(n);
        will_damage = (rand <= p_damage*(es./(IE0+es)));
        if will_damage
            dam(A==new_id)=1;
        end
        % update counter
        ind = ind+1;
    end

    % Biogenesis events
    ids = unique(A(A>0)); % get different ids for mitos
    Es = zeros(1, numel(ids));
    pBiogen = zeros(1, numel(ids));
    for n=1:numel(ids)
        Es(n) = get_stress(ids(n));
        pBiogen(n) = (1-exp(-(lam_bio+Es(n)./(IE0+Es(n)))*dt));
    end

    will_biogen = rand(size(pBiogen))<=pBiogen;
    idf = ids(will_biogen); % ids of those that biogen
    biogen_event(t) =sum(will_biogen); % log biogen events

    ind = 1;
    while numel(idf)>0 && ind<=numel(idf)
        n = idf(ind);
        % remove the index
        idf(idf==n) =[];

        nbrs = get_border(n);
        nbrs_s = get_border_sides(n);
        free_sp = nbrs(A(nbrs)==0);
        free_sp_sides = nbrs_s(A(nbrs_s)==0);

        new_mass = M_single + get_mass(n);

        % if there's no free space or the mass exceed M_max then skip 
        if ~numel(free_sp) || new_mass>M_max
            continue
        end

        % determine the location of the new mito
        if numel(free_sp_sides)>0
            ind_r = randi(numel(free_sp_sides));
            id_bio = free_sp_sides(ind_r);
        else
            ind_r = randi(numel(free_sp));
            id_bio = free_sp(ind_r);
        end

        % Actual biogenesis step
        nn = sum(A(:)==n);
        new_dam = nn*(mean(dam(A==n))+0)/(nn+1);
        A(id_bio)=n;
        dam(A==n) = new_dam;
        dPsi(id_bio)	= dPsi(id_bio)	+ mean(dPsi(A==n));
        H_x(id_bio)		= H_x(id_bio)	+ mean(H_x(A==n));
        NADH_x(id_bio)	= NADH_x(id_bio)+ mean(NADH_x(A==n));
        Pi_x(id_bio)	= Pi_x(id_bio)	+ mean(Pi_x(A==n));
        QH2(id_bio)		= QH2(id_bio)	+ mean(QH2(A==n));
        Cred(id_bio)	= Cred(id_bio)	+ mean(Cred(A==n));
        O2(id_bio)		= O2(id_bio)	+ mean(O2(A==n));
        K_x(id_bio)		= K_x(id_bio)	+ mean(K_x(A==n));
        ATP_x(id_bio)	= ATP_x(id_bio)	+ mean(ATP_x(A==n));
        ADP_x(id_bio)	= ADP_x(id_bio)	+ mean(ADP_x(A==n));
        % update counter
        ind = ind+1;
    end

    % Mitophagy and damage
    %	Only those who die AND have a damage that exceeds a threshold are killed
    dam(dam>=1) = dam(dam>=1)+(1-exp(-k_damage*dt)); % update damage
    ids = unique(A(A>0)); % get different ids for mitos
    which_die = (rand(size(ids))<=p_death); % which die

    high_damaged_mito = ismember(ids, unique(A(dam>=Td)));
    dead = which_die & high_damaged_mito;

    if any(dead)
        dam(ismember(A,ids(dead))) = NaN;
        A(ismember(A,ids(dead))) = 0;
    end
        
    % Solve the PDEs attached to the model
    extra.xN = xN;
    extra.yN = yN;
    extra.A = A;
    extra.dt = dt;
    extra.M_single = M_single;
    extra.M_basal = M_basal;
    extra.XATPase = XATPase;
    
	% Construct diffusion operator
	L = blkdiag(diff_ANP*L_nf, ...	% ATP
			diff_ANP*L_nf, ...	% ADP
			diff_ANP*L_nf, ...	% AMP
			diff_cr*L_nf, ...	% PCr
			diff_cr*L_nf, ...	% Cr
			diff_pi*L_nf, ...	% Pi
			Zh, ...				% ATP_ig
			Zh,	...				% dPsi
			Zh,	...				% H_x
			Zh, ...				% NADH_x
			Zh, ...				% Pi_x
			Zh, ...				% QH2
			Zh, ...				% Cred
			diff_o*L_diri,...	% O2
			Zh, ...				% K_x
			Zh, ...				% ATP_x
			Zh);				% ADP_x

	STATES_INPUT  = [ATP(:); ADP(:); AMP(:); PCr(:); Cr(:); Pi(:); ATP_ig(:); dPsi(:); H_x(:); NADH_x(:); Pi_x(:); QH2(:); Cred(:); O2(:); K_x(:); ATP_x(:); ADP_x(:)];
	pde_out = pdesolve(tmesh(t), STATES_INPUT, L, extra);
   
    % Update state variables
    mitochondria(:,:,t) = A;
    damage(:,:,t)= dam;
    
    % Purify values below machine epsilon to prevent numerical issues
	atp(:,:,t)	  = purify(pde_out.ATP);
	adp(:,:,t)	  = purify(pde_out.ADP);
	amp(:,:,t)	  = purify(pde_out.AMP);
	pcr(:,:,t)	  = purify(pde_out.PCr);
	cr(:,:,t)	  = purify(pde_out.Cr);
	pi(:,:,t)	  = purify(pde_out.Pi);
	atp_ig(:,:,t) = purify(pde_out.ATP_ig);
	dpsi(:,:,t)	  = purify(pde_out.dPsi);
	h_x(:,:,t)	  = purify(pde_out.H_x);
	nadh_x(:,:,t) = purify(pde_out.NADH_x);
	pi_x(:,:,t)	  = purify(pde_out.Pi_x);
	qh2(:,:,t)	  = purify(pde_out.QH2);
	cred(:,:,t)	  = purify(pde_out.Cred);
	o2(:,:,t)	  = purify(pde_out.O2);
	k_x(:,:,t) 	  = purify(pde_out.K_x);
	atp_x(:,:,t)  = purify(pde_out.ATP_x);
	adp_x(:,:,t)  = purify(pde_out.ADP_x);


end
tak1 = min(reshape( mean( mean(o2 , 3), 1) ,1,[])');
tak2 = min(reshape( mean( mean(o2 , 3), 2) ,1,[])');

model_output.tak = [tak1; tak2];

% Store our results
model_output.M 		= mitochondria;
model_output.d 		= damage;
model_output.ATP	= atp;
model_output.ADP	= adp;
model_output.AMP	= amp;
model_output.PCr	= pcr;
model_output.Cr		= cr;
model_output.Pi		= pi;
model_output.ATP_ig	= atp_ig;
model_output.dPsi	= dpsi;
model_output.H_x	= h_x;
model_output.NADH_x	= nadh_x;
model_output.Pi_x	= pi_x;
model_output.QH2	= qh2;
model_output.Cred	= cred;
model_output.O2		= o2;
model_output.K_x	= k_x;
model_output.ATP_x	= atp_x;
model_output.ADP_x	= adp_x;

model_output.fuse 	= fuse_event;
model_output.split 	= split_event;
model_output.biogen = biogen_event;
model_output.stress = stressor_event;

fileName = ['runNumber',num2str(index_input),'_fuse',num2str(par(1)),'_split',num2str(par(2)),'_xatpase',num2str(par(3)),'_o2bdry',num2str(par(4)),'_damage',num2str(par(5)),'.mat'];
save(fileName,'model_output','-v7.3');

% Auxiliary functions
    function id_adj = get_adj(ind_adj)
        adj_list = A(get_border(ind_adj));
        id_adj = (adj_list(adj_list>0));
    end

    function stress = get_stress(id_stress)
        atp_avg = mean(ATP(A==id_stress));
        adp_avg = mean(ADP(A==id_stress));
        atp_ratio = adp_avg./(atp_avg);
        Mx = get_mass(id_stress);
        stress = (atp_ratio./Mx)+get_damage(id_stress);
    end

    function id_nbrs = get_nbrs_indv(id_n)
        [rn,cn]=ind2sub(size(A),id_n);
        id_nbrs =[];
        ret = union(rn+yN*(cn-2),rn+yN*cn);
        ret = union(ret,rn+yN*(cn-1)+1);
        ret = union(ret,rn+yN*(cn-1)-1);
        ret = ret(ret>0 & ret<=xN*yN);
        ret = unique(ret(A(ret)~=id_n));
        if numel(ret)
            id_nbrs = ret';
        end
    end

    function boundary = get_border_sides(id_bd)
        [r_bd,c_bd] = find(A==id_bd);
        boundary =[];
        u = r_bd+yN*(c_bd-2);
        u = union(u,r_bd+yN*(c_bd));
        u = u(u>0 & u<=xN*yN);
        u = u(A(u)~=id_bd);
        u = unique(u);
        if numel(u)
            boundary = u';
        end
    end

    function val = get_new_id(A)
        lb = unique(A(A>0));
        rn = 1:lb(end);
        tmp = rn(~ismember(rn,lb));
        if length(tmp)<1
            val = rn(end)+1;
        else
            val = tmp(1);
        end
    end

    function val = get_damage(id)
        val = mean(dam(A==id));
        if isnan(val)
            val = 0;
        end
    end

    function boundary = get_border(id_bd)
        [r_bd,c_bd] = find(A==id_bd);
        boundary =[];
        ret = r_bd+yN*(c_bd-2);
        ret = union(ret,r_bd+yN*(c_bd));
        ret = union(ret,r_bd+yN*(c_bd-1)+1);
        ret = union(ret,r_bd+yN*(c_bd-1)-1);
        ret = ret(ret>0 & ret<=xN*yN);
        ret = ret(A(ret)~=id_bd);
        ret = unique(ret);
        if numel(ret)
            boundary = ret';
        end
    end

    function mass = get_mass(id_mass)
        mass = M_single*sum(A(:)==id_mass);
    end

    % Real parts and cut of anything below mach_tol to be zero 
    function ret = purify(in)
    	ret = in;
    	ret(~isnan(in)) = real(max(0, in(~isnan(in))-(mach_tol)));
    end
end
