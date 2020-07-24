%Adaption of TestSanWithFunnyCurrentBlockLiteratePaper
function RunToSteadyStateGetCycleLengthAndApd(rExperimentName, pModel,blockLevel,experimentIndex)
    %When solving for a long time with no output we need to increase the default number of steps CVODE can take.
    pModel.MaxSteps= 1e6;
    pModel.t0 = 0;
    pModel.tend = 100*1000;
    pModel.StepSize = 10;%ms
    pModel.Steps = (pMod    
    
    el.tend - pModel.t0 -1)/pModel.StepSize;
    pModel.solution = Solve(pModel.t0,pModel.tend,pModel.Stepsize);% Solve for t=0 to t=100 seconds(100000ms), maximum time step 10ms (will be auto-refined by CVODE).
    %Get output every 0.1ms over a 2 second run for analysing carefully.    
    pModel.t0 = 0;
    pModel.tend = 2*1000;
    pModel.StepSize = 0.1;%ms
    pModel.Steps = Solve(pModel.tend - pModel.t0 -1)/pModel.StepSize;
    pModel.solution = Solve(pModel.t0,pModel.tend,pModel.Stepsize);% Solve for t=0 to t=100 seconds(100000ms), maximum time step 10ms (will be auto-refined by CVODE).

    pModel.voltages = solution.GetVariableAtIndex(voltage_index);
    %CellProperties voltage_properties(voltages, solution.rGetTimes(), -30);
end