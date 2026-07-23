# RSA.jl test suite
#
# IO tests
# Run each tutorial once and check that every input file is read correctly
# Use tutorial 09-HDF5 to test the HDF5 input/output functionality
#

using Test
using Random
using RSA

# Switch in the tutorials directory
tutorial_path = joinpath(@__DIR__,"..","tutorials")
cd(tutorial_path)

@testset "Tutorial 01-Adsorption-Stochastics IO" begin
    inputfile_path = "./01-Adsorption-Stochastics/input.inp"
    Nmolecules, molecules, Ngrids, grids, lattice, events = RSA.read_input(inputfile_path)

    # Molecule input
    @test Nmolecules == 1
    @test molecules[1].label == 1
    @test molecules[1].grids == [1]
    @test molecules[1].fixpointtype == "atoms"
    @test molecules[1].fixpointatoms == [1]
    @test molecules[1].fixpoint ≈ [5.04118; 8.64882; 8.56544;;]
    @test molecules[1].rotationmodus == "angle"
    @test molecules[1].rotationangle == 60.0
    @test molecules[1].rotationvalues == Float64[]
    @test molecules[1].Nrotations == 6
    @test molecules[1].rotations == [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
    @test molecules[1].coordinate_path == "./01-Adsorption-Stochastics/aniline.xyz"
    @test molecules[1].Natoms == 14
    @test molecules[1].dimension == 3
    @test molecules[1].elements == [7, 6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1]
    @test size(molecules[1].coordinates) == (3, 14)
    @test molecules[1].coordinates[:,1] ≈ [5.04118, 8.64882, 8.56544]
    @test length(molecules[1].coordinates_rotated) == 6
    @test size(molecules[1].coordinates_rotated[1]) == (3, 14)
    @test molecules[1].coordinates_rotated[1][:,1] ≈ [0.01747000000000032, -5.28256, 0.5539799999999993]
    @test molecules[1].maxradius ≈ 6.482588887515288
    @test molecules[1].radii == [1.55, 5.39028340002499, 4.118731797306183, 3.0850674079264158, 4.118380393176392, 5.389889747133375, 5.894072637771074, 5.766783945688695, 3.860113147010856, 3.859739115796886, 5.766247817234627, 6.482588887515288, 2.180159486257211, 2.1791210844936386]

    # Grid input
    @test Ngrids == 1
    @test grids[1].label == 1
    @test grids[1].Nuniquepoints == 2
    @test grids[1].uniquepoints ≈ [0.0 1.25942; 0.0 2.18137; 0.0 0.0]
    @test grids[1].dimension == 3
    @test grids[1].Npoints == 1178
    @test size(grids[1].points) == (3, 1178)
    @test grids[1].points[:,1] ≈ [0.0, 0.0, 0.0]
    @test grids[1].points[:,end] ≈ [76.82432, 80.71069, 0.0]
    @test size(grids[1].mapping) == (3, 1178)
    @test grids[1].mapping[:,1] == [1, 0, 0]
    @test grids[1].mapping[:,end] == [2, 30, 18]

    # Lattice input
    @test lattice.vectors ≈ [2.51883 0.0 0.0; 0.0 4.36274 0.0; 0.0 0.0 1.0]
    @test lattice.dimension == 3
    @test lattice.transx == 30
    @test lattice.transy == 18
    @test lattice.Ncellx == 31
    @test lattice.Ncelly == 19
    @test size(lattice.transvectors) == (3, 3)
    @test lattice.transvectors[:,1] ≈ [78.08373, 0.0, 0.0]
    @test size(lattice.inversevectors) == (3, 3)
    @test lattice.inversevectors[:,1] ≈ [0.012806765250584213, 0.0, 0.0]

    # Events input
    @test events.steps == 1000
    @test events.break_steps == true
    @test events.coverage_convergence == 0
    @test events.break_convergence == false
    @test events.Nforce_adsorption == 0
    @test events.force_adsorption == false
    @test events.overlap2d == true
    @test events.overlap3d == false
    @test events.Nevents == 1
    @test events.Nadsorptions == 1
    @test events.Ndiffusions == 0
    @test events.Nrotations == 0
    @test events.Nconformers == 0
    @test events.adsorptions[1].molecule == 1 && events.adsorptions[1].grid == 1 && events.adsorptions[1].weigth == 1.0
end

@testset "Tutorial 03-Multiple-Adsorbates IO" begin
    inputfile_path = "./03-Multiple-Adsorbates/input.inp"
    Nmolecules, molecules, Ngrids, grids, lattice, events = RSA.read_input(inputfile_path)

    # Only new molecule is checked
    @test Nmolecules == 2
    @test molecules[2].label == 2
    @test molecules[2].grids == [1]
    @test molecules[2].fixpointtype == "atoms"
    @test molecules[2].fixpointatoms == [1]
    @test molecules[2].fixpoint ≈ [6.29357; 6.02642; 8.42824;;]
    @test molecules[2].rotationmodus == "angle"
    @test molecules[2].rotationangle == 60.0
    @test molecules[2].Nrotations == 6
    @test molecules[2].rotations == [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
    @test molecules[2].coordinate_path == "./03-Multiple-Adsorbates/pyridine-tilted.xyz"
    @test molecules[2].Natoms == 11
    @test molecules[2].dimension == 3
    @test molecules[2].elements == [7, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1]
    @test size(molecules[2].coordinates) == (3, 11)
    @test molecules[2].coordinates[:,1] ≈ [6.29357, 6.02642, 8.42824]
    @test length(molecules[2].coordinates_rotated) == 6
    @test molecules[2].coordinates_rotated[1][:,1] ≈ [2.1209, -1.8507699999999998, 1.8771499999999985]
    @test molecules[2].maxradius ≈ 4.014883017622579
    @test molecules[2].radii == [1.55, 2.947383170801979, 3.573490753246463, 3.6530742733700627, 3.5699157156674204, 2.9453981445706425, 3.942573332930224, 3.2778629637442407, 4.014883017622579, 4.009211246043985, 3.2772455148344886]
end

@testset "Tutorial 04-Multiple-Grids IO" begin
    inputfile_path = "./04-Multiple-Grids/input.inp"
    Nmolecules, molecules, Ngrids, grids, lattice, events = RSA.read_input(inputfile_path)

    # Only new grid is checked
    @test Ngrids == 2
    @test grids[2].label == 2
    @test grids[2].Nuniquepoints == 4
    @test grids[2].uniquepoints ≈ [0.0 0.0 1.25942 1.25942; 1.45518 2.90908 0.72771 3.63655; 0.0 0.0 0.0 0.0]
    @test grids[2].dimension == 3
    @test grids[2].Npoints == 2356
    @test size(grids[2].points) == (3, 2356)
    @test grids[2].points[:,1] ≈ [0.0, 1.45518, 0.0]
    @test grids[2].points[:,end] ≈ [76.82432, 82.16587, 0.0]
    @test size(grids[2].mapping) == (3, 2356)
    @test grids[2].mapping[:,1] == [1, 0, 0]
    @test grids[2].mapping[:,end] == [4, 30, 18]
end

@testset "Tutorial 05-Rotations IO" begin
    inputfile_path = "./05-Rotations/input.inp"
    Nmolecules, molecules, Ngrids, grids, lattice, events = RSA.read_input(inputfile_path)

    # Only a few checks on the events
    @test events.coverage_convergence == 100
    @test events.break_convergence == true
    @test events.Nforce_adsorption == 25
    @test events.force_adsorption == true
    @test events.Nrotations == 1
    @test events.rotations[1].molecule == 1 && events.rotations[1].grid == 1 && events.rotations[1].weigth == 1.0
end

@testset "Tutorial 06-Diffusions IO" begin
    inputfile_path = "./06-Diffusions/input.inp"
    Nmolecules, molecules, Ngrids, grids, lattice, events = RSA.read_input(inputfile_path)

    # Only a few checks on the events
    @test events.Ndiffusions == 1
    @test events.diffusions[1].molecule == 1
    @test events.diffusions[1].grid_start == 1
    @test events.diffusions[1].grid_end == 1
    @test events.diffusions[1].weigth == 10.0
    @test events.diffusions[1].radius == 2.6
end

@testset "Tutorial 07-Conformer-Changes IO" begin
    inputfile_path = "./07-Conformer-Changes/input.inp"
    Nmolecules, molecules, Ngrids, grids, lattice, events = RSA.read_input(inputfile_path)

    # Only a few checks on the events
    @test events.Nconformers == 1
    @test events.conformers[1].molecule_start == 1
    @test events.conformers[1].molecule_end == 2
    @test events.conformers[1].grid == 1
    @test events.conformers[1].weigth == 10.0
end

@testset "Tutorial 08-Mimic-Reactions IO" begin
    inputfile_path = "./08-Mimic-Reactions/input.inp"
    Nmolecules, molecules, Ngrids, grids, lattice, events = RSA.read_input(inputfile_path)

    # Only a few checks on the events
    @test events.Ndiffusions == 2
    @test events.diffusions[1].molecule == 2 && events.diffusions[1].grid_start == 1 && events.diffusions[1].grid_end == 2 && events.diffusions[1].weigth == 1000.0 && events.diffusions[1].radius == 2.6
    @test events.diffusions[2].molecule == 2 && events.diffusions[2].grid_start == 2 && events.diffusions[2].grid_end == 2 && events.diffusions[2].weigth == 10.0 && events.diffusions[2].radius == 2.6
    @test events.Nconformers == 1
    @test events.conformers[1].molecule_start == 1 && events.conformers[1].molecule_end == 2 && events.conformers[1].grid == 1 && events.conformers[1].weigth == 5.0
end

@testset "Tutorial 09-HDF5 HDF5 write and read test" begin
    Random.seed!(2024)
    inputfile_path = "./09-HDF5/input.inp"
    hdf5_path = "./09-HDF5/input.h5"
    
    # Make sure a leftover file from a previous (failed) run does not block this test
    isfile(hdf5_path) && rm(hdf5_path)

    NRuns = 1
    rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events =
        perform_multiple_rsa_runs(NRuns, inputfile_path; hdf5=true)

    @test isfile(hdf5_path)

    rsa_results2, Nmolecules2, molecules2, Ngrids2, grids2, lattice2, events2 =
        read_hdf5_output_file(hdf5_path)

    # Both sets must be identical
    @test Nmolecules2 == Nmolecules
    @test Ngrids2 == Ngrids
    @test length(rsa_results2) == NRuns
    @test rsa_results2[NRuns].Nsteps == rsa_results[NRuns].Nsteps

    for i in 1:Nmolecules
        @test molecules2[i].label == molecules[i].label
        @test molecules2[i].grids == molecules[i].grids
        @test molecules2[i].Natoms == molecules[i].Natoms
        @test molecules2[i].Nrotations == molecules[i].Nrotations
        @test molecules2[i].coordinates == molecules[i].coordinates
        @test molecules2[i].elements == molecules[i].elements
        @test molecules2[i].maxradius ≈ molecules[i].maxradius
    end

    for i in 1:Ngrids
        @test grids2[i].label == grids[i].label
        @test grids2[i].Npoints == grids[i].Npoints
        @test grids2[i].points == grids[i].points
        @test grids2[i].mapping == grids[i].mapping
    end

    @test lattice2.transx == lattice.transx
    @test lattice2.transy == lattice.transy
    @test lattice2.vectors == lattice.vectors

    @test events2.Nadsorptions == events.Nadsorptions
    @test events2.Nconformers == events.Nconformers

    # Clean up the file this test created
    rm(hdf5_path, force=true)
end