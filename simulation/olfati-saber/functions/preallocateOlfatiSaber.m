function leech = preallocateOlfatiSaber(setup, sim, Tf)
    % Estimate the number of steps based on the final simulation time
    % Preallocate storage for simulation results
    steps = sim.estimateSteps(Tf);
    leech = DataLeech(setup.Agents, steps, 'position', 'velocity', 'u', 'num_N', 'id', 'neighbors');
end