# RSA.jl test suite
#
# Main file
#

using Test
using Random
using RSA

# Vector storing all possible test suits
all_tests = ["rsa", "io"]

# Check whether to run all tests or only a subset (passed as command line arguments)
if isempty(ARGS)
    selected_tests = all_tests
else
	selected_tests = ARGS
end

# Run the selected test suits
@testset verbose=true "RSA test suits" begin
    for ele in selected_tests
        include("$(ele)_tests.jl")
    end
end

