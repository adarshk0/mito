function ret = get_node_value(index,STATES,num_vars,uN)
	ix = index + (0:(num_vars)-1)*uN;
	ret = STATES(ix);
end