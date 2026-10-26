# RSA.jl test suite
#
# PBC tests
# Check the minimum image convention against a brute force search over periodic images
# Check the lattice basis reduction, which allows any choice of transx and transy (e.g. hexagonal cells with transx ≠ transy)
#

using Test
using Random
using RSA

# Switch in the tutorials directory
tutorial_path = joinpath(@__DIR__,"..","tutorials")
cd(tutorial_path)

# Unit cells used below (lattice vectors as columns)
orthogonal_cell = [2.51883 0.0 0.0; 0.0 4.36274 0.0; 0.0 0.0 1.0]
hexagonal_cell_60 = [2.51883 1.25942 0.0; 0.0 2.18137 0.0; 0.0 0.0 1.0]
hexagonal_cell_120 = [2.51883 -1.25942 0.0; 0.0 2.18137 0.0; 0.0 0.0 1.0]
oblique_cell = [3.0 1.1 0.0; 0.0 2.7 0.0; 0.0 0.0 1.0]

# Full lattices: (label, unit cell, transx, transy, basis already reduced?)
test_lattices = [
    ("orthogonal 40x24",    orthogonal_cell,    40, 24, true),
    ("hexagonal 60° 33x33", hexagonal_cell_60,  33, 33, true),
    ("hexagonal 60° 40x24", hexagonal_cell_60,  40, 24, false),
    ("hexagonal 120° 6x1",  hexagonal_cell_120,  6,  1, false),
    ("hexagonal 60° 49x2",  hexagonal_cell_60,  49,  2, false),
    ("oblique 9x4",         oblique_cell,        9,  4, false),
]

# Small helpers (the test environment does not depend on LinearAlgebra)
vector_norm(vector) = sqrt(sum(abs2, vector))
identity_matrix(dimension) = [Int64(i == j) for i in 1:dimension, j in 1:dimension]
inplane_determinant(matrix) = matrix[1,1] * matrix[2,2] - matrix[1,2] * matrix[2,1]

# Shortest periodic image of a cartesian vector by brute force over many images of the full lattice
function brute_force_minimum_image(vector, fulllattice; range = 12)
    shortest = vector
    for shift_x in -range:range, shift_y in -range:range
        candidate = vector + shift_x * fulllattice[:,1] + shift_y * fulllattice[:,2]
        if vector_norm(candidate) < vector_norm(shortest)
            shortest = candidate
        end
    end
    return shortest
end

# Random distance vectors spanning several full lattices
function random_distance_vectors(fulllattice, Nvectors)
    fractional = 6.0 .* (rand(3, Nvectors) .- 0.5)
    fractional[3,:] .= rand(Nvectors) .- 0.5
    return fulllattice * fractional
end

@testset "Lattice basis reduction" begin
    for (label, cell, transx, transy, reduced) in test_lattices
        fulllattice, invfulllattice = RSA.replicate_lattice_with_translation(cell, transx, transy)
        reducedlattice, invreducedlattice, reductionmatrix = RSA.reduce_lattice_basis(fulllattice, invfulllattice)

        @test RSA.is_reduced_lattice_basis(fulllattice) == reduced
        @test RSA.is_reduced_lattice_basis(reducedlattice)

        # A reduced basis is returned unchanged
        if reduced
            @test reducedlattice == fulllattice
            @test invreducedlattice == invfulllattice
            @test reductionmatrix == identity_matrix(3)
        end

        # The reduced basis spans the same periodic images: unimodular integer transformation
        @test abs(inplane_determinant(reductionmatrix)) == 1
        @test reductionmatrix[:,3] == [0, 0, 1] && reductionmatrix[3,:] == [0, 0, 1]
        @test reducedlattice * reductionmatrix ≈ fulllattice
        @test reducedlattice[:,3] == fulllattice[:,3]
    end
end

@testset "Minimum image convention vs. brute force" begin
    Random.seed!(2026)
    for (label, cell, transx, transy, reduced) in test_lattices
        fulllattice, invfulllattice = RSA.replicate_lattice_with_translation(cell, transx, transy)
        reducedlattice, invreducedlattice, reductionmatrix = RSA.reduce_lattice_basis(fulllattice, invfulllattice)
        distance_vectors = random_distance_vectors(fulllattice, 2000)

        # Matrix version (cartesian input)
        mic_vectors = RSA.apply_minimum_image_convention(distance_vectors, reducedlattice, invreducedlattice)

        Nshortest = 0
        Nperiodic = 0
        Nfractional = 0
        for column_id in axes(distance_vectors, 2)
            vector = distance_vectors[:, column_id]
            reference = brute_force_minimum_image(vector, fulllattice)

            # Same length as the shortest image
            Nshortest += isapprox(vector_norm(mic_vectors[:, column_id]), vector_norm(reference); atol = 1e-8)

            # Result is a periodic image of the input (differs by a lattice vector of the full lattice)
            shift = invfulllattice * (mic_vectors[:, column_id] - vector)
            Nperiodic += all(abs.(shift .- round.(shift)) .< 1e-8)

            # Fractional input gives the same result
            fractional_result = RSA.minimum_image_from_fractional_coordinates(invfulllattice * vector, reducedlattice, reductionmatrix)
            Nfractional += isapprox(fractional_result, mic_vectors[:, column_id]; atol = 1e-8)
        end
        @test Nshortest == size(distance_vectors, 2)
        @test Nperiodic == size(distance_vectors, 2)
        @test Nfractional == size(distance_vectors, 2)
    end
end

@testset "Neighbour images of a non-reduced basis are not sufficient" begin
    # Documents why the reduction is necessary: for a strongly elongated full lattice (here a hexagonal
    # 50x3 cell), searching only the neighbouring images of the original basis misses the shortest image
    # for many distance vectors. Mildly skewed cells (e.g. hexagonal 41x25) happen to work without the
    # reduction, but only the reduced basis guarantees the shortest image for every cell.
    Random.seed!(2026)
    fulllattice, invfulllattice = RSA.replicate_lattice_with_translation(hexagonal_cell_60, 49, 2)
    distance_vectors = random_distance_vectors(fulllattice, 2000)
    Nmissed = 0
    for column_id in axes(distance_vectors, 2)
        vector = distance_vectors[:, column_id]
        unreduced_result = RSA.minimum_image_vector(vector, fulllattice, invfulllattice)
        Nmissed += vector_norm(unreduced_result) > vector_norm(brute_force_minimum_image(vector, fulllattice)) + 1e-8
    end
    @test Nmissed > 0
end

@testset "Tutorial 11-Hexagonal-Cell with transx ≠ transy" begin
    # The input reader generates the reduced basis
    Nmolecules, molecules, Ngrids, grids, lattice, events = RSA.read_input("./11-Hexagonal-Cell/input.inp")
    @test lattice.reducedcellvectors == lattice.transcellvectors
    @test lattice.reductionmatrix == identity_matrix(3)

    # Same input with unequal translations (written next to the tutorial input, removed afterwards)
    inputfile_path = "./11-Hexagonal-Cell/input-unequal-translations.inp"
    input_text = read("./11-Hexagonal-Cell/input.inp", String)
    input_text = replace(input_text, "transx = 33" => "transx = 40", "transy = 33" => "transy = 24")
    write(inputfile_path, input_text)

    Random.seed!(2024)
    NRuns = 20
    rsa_results_unequal, Nmolecules, molecules, Ngrids, grids, lattice, events = perform_multiple_rsa_runs(NRuns, inputfile_path)
    @test lattice.transx == 40
    @test lattice.transy == 24
    @test lattice.reductionmatrix != identity_matrix(3)

    Random.seed!(2024)
    rsa_results_equal, = perform_multiple_rsa_runs(NRuns, "./11-Hexagonal-Cell/input.inp")

    # The number of (point, rotation) combinations blocked by the first adsorbate only depends on the
    # molecule and the local lattice, not on the shape of the full lattice
    blocked_by_first_adsorbate(result) = result.stepinfo[1,1] - result.stepinfo[1,2]
    reference = blocked_by_first_adsorbate(rsa_results_equal[1])
    @test all(blocked_by_first_adsorbate(result) == reference for result in rsa_results_equal)
    @test all(blocked_by_first_adsorbate(result) == reference for result in rsa_results_unequal)
    @test rsa_results_unequal[1].stepinfo[1,1] == 41 * 25 * 6

    # Clean up the file this test created
    rm(inputfile_path, force=true)

end
