function u = StrangSplitting(t0, u0, lin, arg)
    % Load up other stuff
    dt = arg.dt;
    u  = u0;
    % Non-linear part is solved with ode23s
    % Linear part is solved using ode15s but with 'MaxOrder' being 2
    u = solve_lin([t0, t0+0.5*dt], u, lin); 
    u = solve_nonlin([t0, t0+dt], u, arg);
    u = solve_lin([t0+0.5*dt, t0+dt], u, lin);
end