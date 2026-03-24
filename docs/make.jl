push!(LOAD_PATH,"../src/")
using Documenter
using RSA

molecule_subsection_silanes = ["systems/2026-tripropoxymethylsilane.md"]
molecule_subsection_acids = ["systems/2026-methanesulfonicacid.md"]
molecule_subsection_thiols = ["systems/2026-diethylsulfide.md", "systems/2026-diisopropylsulfide.md", "systems/2026-dipropylsulfide.md"]
molecule_subsections_cyclic = ["systems/2026-aniline.md", "systems/2026-pyridine.md", "systems/2026-pyrrole.md"]
surface_subsection = ["systems/2026-Cu.md", "systems/2026-SiO2.md"]
tutorial_subsection = ["tutorials/01-adsorption-stochastics.md", "tutorials/02-benchmarks.md", "tutorials/03-multiple-adsorbates.md", "tutorials/04-multiple-grids.md", "tutorials/05-events-rotations.md", "tutorials/06-events-diffusions.md", "tutorials/07-events-conformer-changes.md", "tutorials/08-events-conversions.md", "tutorials/09-hdf5.md"]

# Local non-ideal solution
#makedocs(sitename="RSA.jl", remotes=nothing,
   #format = Documenter.HTML(assets = ["assets/custom.css"], prettyurls = false),

makedocs(sitename="RSA.jl",
   format = Documenter.HTML(assets = ["assets/custom.css"], prettyurls = true),
   pages = ["Home" => "index.md",
            "About RSA" => "rsa.md",
            "Installation" => "install.md",
            "Input Keywords" => ["Keyword Types" => "keywords.md",
                                 "Molecules" => "molecules.md",
                                 "Lattice" => "lattice.md",
                                 "Grids" => "grids.md",
                                 "Events" => "events.md"],
            "Running Simulations" => "run.md",
            hide("Tutorials" => "tutorials/tutorials.md", tutorial_subsection),
            "Systems" => [hide("Molecule Library" => "systems/molecules.md", [molecule_subsection_silanes; molecule_subsection_thiols; molecule_subsections_cyclic; molecule_subsection_acids]),
                          hide("Surface Library" => "systems/surfaces.md", surface_subsection)],
            "Developers" => "developers.md"
   ]
)

# comment for local build
deploydocs(
    repo = "github.com/Tonner-Zech-Group/RSA.jl.git",
)

