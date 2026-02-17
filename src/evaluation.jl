# Functions for plotting and evaluation of data

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


function check_kMC_convergence(final_occupations::Array{Int64}, times::Vector{Float64}, tolerance::Float64)
    
    # Read the last elemet of every state and calculate the tolerance
    array_size = size(final_occupations,1)
    last_element = final_occupations[array_size,:]
    Nstates = size(final_occupations,2)
    tolerance_range = Array{Float64}(undef,2,Nstates)
    for i in eachindex(last_element)
        tolerance_range[1,i] = last_element[i] - last_element[i] * tolerance
        tolerance_range[2,i] = last_element[i] + last_element[i] * tolerance
    end

    # Check when the tolerance is reached
    convergence_element = 0
    convergence_time = 0.0
    convergence_reached = false
    for kmc_step in axes(final_occupations, 1)

        for state in axes(final_occupations, 2)

            if tolerance_range[1,state] <= final_occupations[kmc_step,state] <= tolerance_range[2,state]
                convergence_reached = true
            else
                convergence_reached = false
                break
            end

            if convergence_reached == true
                convergence_element = kmc_step
                convergence_time = times[kmc_step]
            end

        end

        if convergence_reached == true
            break
        end

    end

    # Return the result
    return convergence_reached, convergence_element, convergence_time, tolerance_range

end