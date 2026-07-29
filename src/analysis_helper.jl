#
# General head section of any file within this module 
#

#
# Using statements
#

#
# Include statements
#

#
# Definition of global variables
#

#
# Specific section for this file
#

# A helper function to define the y-axes label of histograms
function histogram_label(normalized)

    if normalized == true
        return "Normalized count"
    else
        return "Count"
    end

end

# A function to calculate the scaling factor for gaussian distributions so that they fit to a histogram
# The area of the visible bins (count x width) has to fit the gaussian distribution (usually normalized to 1)
function gaussian_scaling_factor(data_vector, bins)

    # Count all values within the range of the bins
    Nbinned = 0
    for value in data_vector
        if value ≥ first(bins) && value ≤ last(bins)
            Nbinned += 1
        end
    end

    # Return the scaling factor: number of shown values times the bin width
    return Float64(Nbinned) * Float64(step(bins))

end