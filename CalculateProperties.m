function CalculateProperties(mrTime, mrVoltage,DBL_MAX,prev_t, mThreshold)
    if length(mrTime)<1
        message("Insufficient time steps to calculate physiological properties.")
    end
    if length(mrTime)!=length(mrVoltage)
        message("Time and Voltage series should be the same length.")
    end
    max_upstroke_velocity = -DBL_MAX;
    current_time_of_upstroke_velocity = 0;
    current_resting_value = DBL_MAX;
    current_minimum_velocity = DBL_MAX;
    current_peak_time = -DBL_MAX;
    prev_voltage_derivative = 0;
    mRestingValues = mrVoltage;
    ap_counter = 0;
    counter_of_plateau_depolarisations = 0;
    switching_phase = false;
    found_a_flat_bit = false;
    ap_phase = 'BELOWTHRESHOLD';
    
    time_steps = length(mrTime) - 1;
    v = mrVoltage(0);
    t = mrTime(0);
    prev_v =v;
 %   voltage_derivative
    resting_potential_gradient_threshold = 1e-2;
    
    for i=1:time_steps
        v = mrVoltage(i);
        t = mrTime(i);
        if t==prev_t
            voltage_derivative = 0.0;
        else
            voltage_derivative = (v - prev_v) / (t - prev_t);
        end
    
        % Look for the max upstroke velocity and when it happens (could be below or above threshold).
        if voltage_derivative >= max_upstroke_velocity
            max_upstroke_velocity = voltage_derivative;
            current_time_of_upstroke_velocity = t;
        end
        
        switch (ap_phase)
            case BELOWTHRESHOLD
                if abs(voltage_derivative) <= current_minimum_velocity & abs(voltage_derivative) <= resting_potential_gradient_threshold)
                    current_minimum_velocity = abs(voltage_derivative);
                    current_resting_value = prev_v;
                    found_a_flat_bit = true;
                elseif prev_v < current_resting_value & found_a_flat_bit==false
                    current_resting_value = prev_v;
                end
                % If cross the threshold, this counts as an AP
                if v > mThreshold & prev_v <= mThreshold
                    % register the resting value and re-initialise the minimum velocity
                    mRestingValues(i)= current_resting_value;                
                    current_minimum_velocity = DBL_MAX;
                    current_resting_value = DBL_MAX;
                    
                    % Register the onset time. Linear interpolation.
                    mOnsets(i) = (prev_t+(t-prev_t)/(v-prev_v)*(mThreshold - prev_v));
                    %If it is not the first AP, calculate cycle length for the last two APs
                    if ap_counter > 0
                         mCycleLengths(i) = (mOnsets(ap_counter) - mOnsets(ap_counter - 1));
                    end
                    % Re-initialise max_upstroke_velocity
                    max_upstroke_velocity = -DBL_MAX;
                    
                    switching_phase = true;
                    found_a_flat_bit = false;
                    ap_phase = ABOVETHRESHOLD;
                else
                    break;
                end
                %no break here - deliberate fall through to next case
                
            case ABOVETHRESHOLD
                if v > current_peak
                    current_peak = v;
                    current_peak_time = t;
                end
                
                %we check whether we have above threshold depolarisation
                %and only if if we haven't just switched from below threshold at this time step
                %The latter is to avoid recording things depending on resting behaviour (in case of sudden upstroke from rest)
                if prev_voltage_derivative <= 0 & voltage_derivative > 0 & switching_phase != false    
                    counter_of_plateau_depolarisations = counter_of_plateau_depolarisations+1;
                end
                % From the next time step, we are not "switching phase" any longer
                % we want to check for above threshold deolarisations
                switching_phase = false;
                if v < mThreshold & prev_v >= mThreshold
                    %Register peak value for this AP
                    mPeakValues(i) = current_peak;
                    mTimesAtPeakValues(i) = current_peak_time;    
                    %Re-initialise the current_peak
                    current_peak = mThreshold;
                    current_peak_time = -DBL_MAX;
                    
                    %Register maximum upstroke velocity for this AP
                    mMaxUpstrokeVelocities(i) = max_upstroke_velocity;
                    %Re-initialise max_upstroke_velocity
                    max_upstroke_velocity = -DBL_MAX;
  
                    
                    %Register time when maximum upstroke velocity occurred for this AP
                    mTimesAtMaxUpstrokeVelocity(i) = current_time_of_upstroke_velocity;   
                    % Re-initialise current_time_of_upstroke_velocity=t;        
                    current_time_of_upstroke_velocity = 0.0;
                    
                    mCounterOfPlateauDepolarisations(i) = counter_of_plateau_depolarisations;
                    
                    %update the counters
                    ap_counter = ap_counter+1;
                    ap_phase = BELOWTHRESHOLD;
                    
                    %reinitialise counter of plateau depolarisations
                    counter_of_plateau_depolarisations = 0;
                end        
        end
        prev_v = v;
        prev_t = t;
        prev_voltage_derivative = voltage_derivatage;
    end
    
    % One last check. If the simulation ends halfway through an AP
    % i.e. if the vectors of onsets has more elements than the vectors
    % of peak and upstroke properties (that are updated at the end of the AP),
    % then we register the peak and upstroke values so far
    % for the last incomplete AP.
    if length(mOnsets) > length(mMaxUpstrokeVelocities)
        mMaxUpstrokeVelocities(i) =max_upstroke_velocity;        
        mPeakValues(i) = current_peak;
        mTimesAtPeakValues(i) = current_peak_time;
        mTimesAtMaxUpstrokeVelocity(i) = current_time_of_upstroke_velocity;
        mUnfinishedActionPotentials = true; 
    end   
end

