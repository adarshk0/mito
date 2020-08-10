function yout = solve_nonlin(tspan, STATES, extra)
	A = extra.A;
	xN = extra.xN;
	yN = extra.yN;
    M_single = extra.M_single;
    M_basal = extra.M_basal;
    XATPase = extra.XATPase;
    
	id_myo = find(A==0);
	id_mit = find(A~=0);

	num_myo  = numel(id_myo);
	num_mit = numel(id_mit);

	id_myo_to_solve = sparse(num_myo, num_myo);
	id_mit_to_solve = sparse(num_mit, num_mit);

	tol = 1e-9;

	uN = yN*xN;
	num_vars = numel(STATES)/uN;
	yin = reshape(STATES,uN,num_vars);
	yout = sparse(uN,num_vars);
    
    % Find ids of mitos that are on the boundary 
    bc=sparse(yN,xN);
    bc(1,:) = 1; bc(end,:) = 1; 
    bc(:,1)=1; bc(:,end) =1;
    ix_bdry = find(bc & A);
    
	% Find id's of myo elements that are below tol
	tmp = id_myo;
	n = 1;
	while numel(tmp)>0
		ix  = tmp(1);
		ics = get_node_value(ix,STATES,num_vars,uN)';
		% Calculate other nodes with the same inital value
		also_solved = find(vecnorm(yin-ics,2,2)<tol);

		% Merge ids
      	ix = union(ix, also_solved);
		tmp = setdiff(tmp, ix);

		% Save the ids that are similar
		id_myo_to_solve(n,:) = sparse(1,1:numel(ix),ix,1,num_myo);
		n = n+1;
	end

	% Find the mass of all mitos
	mx_mass = zeros(1,num_mit);
	parfor n=1:num_mit
        id = A(id_mit(n));
		mx_mass(n) = M_single*sum(A(:)==id);
    end
    
    tmp = id_mit;
    tmp = setdiff(tmp, find(ix_bdry)); % Remove boundary mitos
	n = 1;
	while numel(tmp)>0
		ix  = tmp(1);
		mass_ix = mx_mass(find(id_mit==ix,1));
		ids_same_mass = find(mx_mass==mass_ix);
		ics = get_node_value(ix, STATES,num_vars,uN)';
		% Calculate other nodes with the same inital value
        
		also_solved = find(vecnorm(yin-ics,2,2)<tol);
		also_solved = intersect(also_solved,ids_same_mass);
		% Merge ids
		ix = union(ix, also_solved);
		tmp = setdiff(tmp, ix);
		% Save the ids that are similar
		id_mit_to_solve(n,:) = sparse(1,1:numel(ix),ix,1,num_mit);
		n = n+1;
	end

	nnz_mit = nnz(id_mit_to_solve(:,1));
	nnz_myo = nnz(id_myo_to_solve(:,1));
    num_bdry_mit = numel(ix_bdry);
    
    y_mit = zeros(nnz_mit, num_vars);
    y_mit_bdry = zeros(num_bdry_mit, num_vars);
    y_myo = zeros(nnz_myo, num_vars);
    
    % Solve mitochondria on the boundary 
	parfor m=1:num_bdry_mit
		ix = id_mit_to_solve(m,:);
		ix = full(ix(ix~=0));
		mass = mx_mass(find(id_mit==ix(1),1));
		ics = get_node_value(ix(1),STATES,num_vars,uN);
		y_mit_bdry(m,:) = main_mito(tspan, ics, M_basal, mass, true);
    end
    
    % Solve mitochondria that are NOT on the boundary
	parfor m=1:nnz_mit
		ix = id_mit_to_solve(m,:);
		ix = full(ix(ix~=0));
		mass = mx_mass(find(id_mit==ix(1),1));
		ics = get_node_value(ix(1),STATES,num_vars,uN);
		y_mit(m,:) = main_mito(tspan, ics, M_basal, mass, false);
    end
    
    % Solve myofibril equations
	parfor m=1:nnz_myo
		ix = id_myo_to_solve(m,:);
		ix = full(ix(ix~=0));
		ics = get_node_value(ix(1),STATES,num_vars,uN);
		y_myo(m,:) = main_myo(tspan, ics, XATPase);
    end
    
    for m=1:num_bdry_mit
        ix = ix_bdry(m);
        yout(ix,:) = y_mit_bdry(m,:).*ones(numel(ix),num_vars);
    end
    
    for m=1:nnz_mit
        ix = id_mit_to_solve(m,:);
		ix = full(ix(ix~=0));
        yout(ix,:) = y_mit(m,:).*ones(numel(ix),num_vars);
    end
    
    for m=1:nnz_myo
        ix = id_myo_to_solve(m,:);
		ix = full(ix(ix~=0));
        yout(ix,:) = y_myo(m,:).*ones(numel(ix),num_vars);
    end
end