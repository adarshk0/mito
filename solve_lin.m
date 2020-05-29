function u = solve_lin(tspan, ics, lin)
    opts = odeset('Jacobian', lin);
    [~,u]=ode15s(@(t,y) lin*y, tspan, ics, opts);
    u = u(end,:).';
end