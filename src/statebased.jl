# Collection of functions focusing on state based kMC

#
# General head section of any file within this module 
#

#
# Include statements
#

#
# Export statements
#
 
#
# Definition of global variables
#

#
# Specific section for this file
#

@doc raw"""

    time_progression(total_time::Float64, total_rate_constant::Float64)

Based on a random number `r` (in the range 0 to 1) and the total rate constant `k` the time `t` progression is calculated based on:
```math
t_{old} = t_{new} - \frac{ln(1 - r)}{k}
```
Returning the value of ``t_{old}``.
"""
function time_progression(total_time::Float64, total_rate_constant::Float64)
    random_number = 1.0 - rand()
    delta_time = -log(random_number)/total_rate_constant
    total_time = total_time + delta_time
    return total_time
end


@doc raw"""

    create_event_list!(occupation,event_list,event_possible)

This is the slowest implementation to create a list of possible events. Here, the `occupation` vector is used to find occupied states.
This information is used to set all events in the `event_possible` vector to `true`` if the `event_list` vector is indicating that the event is starting from 
the occupied state.
"""
function create_event_list!(occupation,event_list,event_possible)
    
    # Check whether a state is occupied
    for i in 1:length(occupation)
        if occupation[i] > 0

            # Check for all events starting from that state
            for j in 1:size(event_list,1)
                if event_list[j,1] == i
                    event_possible[j] = true
                end
            end
        
        else

            # Check for all events starting from that state
            for j in 1:size(event_list,1)
                if event_list[j,1] == i
                    event_possible[j] = false
                end
            end
            
        end
    end
end