function sweep_timestep()

    for timeStep = 0.01:0.01:0.05
        skeleton_cpsclass_reach(2, timeStep);
    end

end