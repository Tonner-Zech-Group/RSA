# RSA.jl test suite
#
# Full RSA runs
# Run each tutorial once and check that the results are consistent with previous runs
#

using Test
using Random
using RSA

# Switch in the tutorials directory
tutorial_path = joinpath(@__DIR__,"..","tutorials")
cd(tutorial_path)

@testset "Tutorial 01-Adsorption-Stochastics" begin
    Random.seed!(2024)
    inputfile_path = "./01-Adsorption-Stochastics/input.inp"
    NRuns = 1000
    rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events = perform_multiple_rsa_runs(NRuns, inputfile_path);
    steps = rsa_results[1000].Nsteps

    @test rsa_results[1000].Nsteps == 70
    @test rsa_results[1000].Nevents == reshape([[69, 0, 0, 0]], 1, 1)
    @test rsa_results[1000].status[:,1] == [1, 1, 360, 6]
    @test rsa_results[1000].status[:,steps - 1] == [1, 1, 28, 2]
    @test rsa_results[1000].randomseed[1] ≈ 0.3055841475195493
    @test rsa_results[1000].randomseed[steps - 1] ≈ 0.6677994071843252
    @test rsa_results[1000].stepinfo[:,1] == [7068, 7068, 0, 0, 0, 1, 360, 1, 1, 0, 6, 0]
    @test rsa_results[1000].stepinfo[:,steps - 1] == [1, 1, 0, 0, 0, 1, 28, 1, 1, 0, 2, 0]
end

@testset "Tutorial 02-Benchmarks CellSize 20" begin
    Random.seed!(2024)
    inputfile_path = "./02-Benchmarks/CellSize/input-20.inp"
    NRuns = 1000
    rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events = perform_multiple_rsa_runs(NRuns, inputfile_path);
    steps = rsa_results[1000].Nsteps

    @test rsa_results[1000].Nsteps == 12
    @test rsa_results[1000].Nevents == reshape([[11, 0, 0, 0]], 1, 1)
    @test rsa_results[1000].status[:,1] == [1, 1, 60, 6]
    @test rsa_results[1000].status[:,steps - 1] == [1, 1, 64, 1]
    @test rsa_results[1000].randomseed[1] ≈ 0.3888317889849616
    @test rsa_results[1000].randomseed[steps - 1] ≈ 0.35161754638131215
    @test rsa_results[1000].stepinfo[:,1] == [924, 924, 0, 0, 0, 1, 60, 1, 1, 0, 6, 0]
    @test rsa_results[1000].stepinfo[:,steps - 1] == [1, 1, 0, 0, 0, 1, 64, 1, 1, 0, 1, 0]
end

@testset "Tutorial 02-Benchmarks CellSize 80" begin
    Random.seed!(2024)
    inputfile_path = "./02-Benchmarks/CellSize/input-80.inp"
    NRuns = 1
    rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events = perform_multiple_rsa_runs(NRuns, inputfile_path);
    steps = rsa_results[1].Nsteps

    @test rsa_results[1].Nsteps == 124
    @test rsa_results[1].Nevents == reshape([[123, 0, 0, 0]], 1, 1)
    @test rsa_results[1].status[:,1] == [1, 1, 1827, 6]
    @test rsa_results[1].status[:,steps - 1] == [1, 1, 1963, 5]
    @test rsa_results[1].randomseed[1] ≈ 0.8911575406016274
    @test rsa_results[1].randomseed[steps - 1] ≈ 0.7080784890612472
    @test rsa_results[1].stepinfo[:,1] == [12300, 12300, 0, 0, 0, 1, 1827, 1, 1, 0, 6, 0]
    @test rsa_results[1].stepinfo[:,steps - 1] == [5, 5, 0, 0, 0, 1, 1963, 1, 1, 0, 5, 0]
end

@testset "Tutorial 02-Benchmarks Rotations 5" begin
    Random.seed!(2024)
    inputfile_path = "./02-Benchmarks/Rotations/input-5.inp"
    NRuns = 1
    rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events = perform_multiple_rsa_runs(NRuns, inputfile_path);
    steps = rsa_results[1].Nsteps

    @test rsa_results[1].Nsteps == 84
    @test rsa_results[1].Nevents == reshape([[83, 0, 0, 0]], 1, 1)
    @test rsa_results[1].status[:,1] == [1, 1, 1050, 57]
    @test rsa_results[1].status[:,steps - 1] == [1, 1, 685, 35]
    @test rsa_results[1].randomseed[1] ≈ 0.8911575406016274
    @test rsa_results[1].randomseed[steps - 1] ≈ 0.5261756712458538
    @test rsa_results[1].stepinfo[:,1] == [84816, 84816, 0, 0, 0, 1, 1050, 1, 1, 0, 57, 0]
    @test rsa_results[1].stepinfo[:,steps - 1] == [1, 1, 0, 0, 0, 1, 685, 1, 1, 0, 35, 0]
end

@testset "Tutorial 02-Benchmarks Rotations 360" begin
    Random.seed!(2024)
    inputfile_path = "./02-Benchmarks/Rotations/input-360.inp"
    NRuns = 1000
    rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events = perform_multiple_rsa_runs(NRuns, inputfile_path);
    steps = rsa_results[1000].Nsteps

    @test rsa_results[1000].Nsteps == 65
    @test rsa_results[1000].Nevents == reshape([[64, 0, 0, 0]], 1, 1)
    @test rsa_results[1000].status[:,1] == [1, 1, 201, 1]
    @test rsa_results[1000].status[:,steps - 1] == [1, 1, 563, 1]
    @test rsa_results[1000].randomseed[1] ≈ 0.17018977664335366
    @test rsa_results[1000].randomseed[steps - 1] ≈ 0.20809832180717325
    @test rsa_results[1000].stepinfo[:,1] == [1178, 1178, 0, 0, 0, 1, 201, 1, 1, 0, 1, 0]
    @test rsa_results[1000].stepinfo[:,steps - 1] == [1, 1, 0, 0, 0, 1, 563, 1, 1, 0, 1, 0]
end

@testset "Tutorial 03-Multiple-Adsorbates" begin
    Random.seed!(2024)
    inputfile_path = "./03-Multiple-Adsorbates/input.inp"
    NRuns = 1
    rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events = perform_multiple_rsa_runs(NRuns, inputfile_path);
    steps = rsa_results[1].Nsteps

    @test rsa_results[1].Nsteps == 153
    @test rsa_results[1].Nevents == reshape([[98, 0, 0, 0, ],[54, 0, 0, 0]], 2, 1)
    @test rsa_results[1].status[:,1] == [2, 1, 1050, 4]
    @test rsa_results[1].status[:,steps - 1] == [1, 1, 400, 1]
    @test rsa_results[1].randomseed[1] ≈ 0.8911575406016274
    @test rsa_results[1].randomseed[steps - 1] ≈ 0.10219795741731885
    @test rsa_results[1].stepinfo[:,1] == [14136, 14136, 0, 0, 0, 1, 1050, 2, 1, 0, 4, 0]
    @test rsa_results[1].stepinfo[:,steps - 1] == [2, 2, 0, 0, 0, 1, 400, 1, 1, 0, 1, 0]
end

@testset "Tutorial 04-Multiple-Grids" begin
    Random.seed!(2024)
    inputfile_path = "./04-Multiple-Grids/input.inp"
    NRuns = 1
    rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events = perform_multiple_rsa_runs(NRuns, inputfile_path);
    steps = rsa_results[1].Nsteps

    @test rsa_results[1].Nsteps == 122
    @test rsa_results[1].Nevents == reshape([[0, 0, 0, 0], [39, 0, 0, 0], [82, 0, 0, 0], [0, 0, 0, 0]], 2, 2)
    @test rsa_results[1].status[:,1] == [1, 2, 1844, 1]
    @test rsa_results[1].status[:,steps - 1] == [1, 2, 741, 1]
    @test rsa_results[1].randomseed[1] ≈ 0.8911575406016274
    @test rsa_results[1].randomseed[steps - 1] ≈ 0.18096461991769086
    @test rsa_results[1].stepinfo[:,1] == [14136, 14136, 0, 0, 0, 2, 1844, 1, 1, 0, 1, 0]
    @test rsa_results[1].stepinfo[:,steps - 1] == [3, 3, 0, 0, 0, 2, 741, 1, 1, 0, 1, 0]
end

@testset "Tutorial 05-Rotations" begin
    Random.seed!(2024)
    inputfile_path = "./05-Rotations/input.inp"
    NRuns = 1
    rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events = perform_multiple_rsa_runs(NRuns, inputfile_path);
    steps = rsa_results[1].Nsteps

    @test rsa_results[1].Nsteps == 432
    @test rsa_results[1].Nevents == reshape([[80, 352, 0, 0]], 1, 1)
    @test rsa_results[1].status[:,1] == [1, 1, 1050, 2]
    @test rsa_results[1].status[:,79] == [1, 1, 382, 3]
    @test rsa_results[1].randomseed[1] ≈ 0.8911575406016274
    @test rsa_results[1].randomseed[steps - 1] ≈ 0.2722254094134665
    @test rsa_results[1].stepinfo[:,1] == [7068, 7068, 0, 0, 0, 1, 1050, 1, 1, 0, 5, 0]
    @test rsa_results[1].stepinfo[:,steps - 1] == [29, 0, 29, 0, 0, 1, 257, 1, 2, 0, 1, 0]
end

@testset "Tutorial 06-Diffusions" begin
    Random.seed!(2024)
    inputfile_path = "./06-Diffusions/input.inp"
    NRuns = 1
    rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events = perform_multiple_rsa_runs(NRuns, inputfile_path);
    steps = rsa_results[1].Nsteps

    @test rsa_results[1].Nsteps == 1299
    @test rsa_results[1].Nevents == reshape([[143, 0, 1156, 0]], 1, 1)
    @test rsa_results[1].status[:,1] == [1, 1, 2021, 2]
    @test rsa_results[1].status[:,142] == [1, 1, 1429, 2]
    @test rsa_results[1].randomseed[1] ≈ 0.8911575406016274
    @test rsa_results[1].randomseed[steps - 1] ≈ 0.7244166611794907
    @test rsa_results[1].stepinfo[:,1] == [7068, 7068, 0, 0, 0, 1, 2100, 1, 1, 0, 2, 0]
    @test rsa_results[1].stepinfo[:,steps - 1] == [173, 0, 0, 173, 0, 1, 1746, 1, 3, 1, 1672, 2]
end

@testset "Tutorial 07-Conformer-Changes" begin
    Random.seed!(2024)
    inputfile_path = "./07-Conformer-Changes/input.inp"
    NRuns = 1
    rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events = perform_multiple_rsa_runs(NRuns, inputfile_path);
    steps = rsa_results[1].Nsteps

    @test rsa_results[1].Nsteps == 179
    @test rsa_results[1].Nevents == reshape([[81, 0, 0, 47], [50, 0, 0, 0]], 2, 1)
    @test rsa_results[1].status[:,1] == [2, 1, 1050, 4]
    @test rsa_results[1].status[:,130] == [1, 1, 975, 3]
    @test rsa_results[1].randomseed[1] ≈ 0.8911575406016274
    @test rsa_results[1].randomseed[steps - 1] ≈ 0.19279484443680805
    @test rsa_results[1].stepinfo[:,1] == [14136, 14136, 0, 0, 0, 1, 1050, 2, 1, 0, 4, 0]
    @test rsa_results[1].stepinfo[:,steps - 1] == [2, 2, 0, 0, 0, 1, 213, 1, 1, 0, 3, 0]
end

@testset "Tutorial 08-Mimic-Reactions" begin
    Random.seed!(2024)
    inputfile_path = "./08-Mimic-Reactions/input.inp"
    NRuns = 1
    rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events = perform_multiple_rsa_runs(NRuns, inputfile_path);
    steps = rsa_results[1].Nsteps

    @test rsa_results[1].Nsteps == 2500
    @test rsa_results[1].Nevents == reshape([[186, 0, 0, 183],  [0, 0, 183, 0], [0, 0, 0, 0], [0, 0, 1948, 0]], 2, 2)
    @test rsa_results[1].status[:,1] == [2, 2, 2100, 3]
    @test rsa_results[1].status[:,185] == [1, 1, 2844, 4]
    @test rsa_results[1].randomseed[1] ≈ 0.8911575406016274
    @test rsa_results[1].randomseed[steps - 1] ≈ 0.8334743116905989
    @test rsa_results[1].stepinfo[:,1] == [56544, 56544, 0, 0, 0, 1, 3150, 1, 1, 0, 3, 0]
    @test rsa_results[1].stepinfo[:,steps - 1] == [286, 0, 0, 281, 5, 2, 2005, 2, 3, 2, 1929, 1]
end