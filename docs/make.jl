using Documenter, Latlib

makedocs(
    sitename = "Latlib.jl",
    modules = [Latlib],
    authors = "Alexander Wietek and contributors",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://awietek.github.io/Latlib.jl",
        edit_link = "main",
    ),
    pages = [
        "Home" => "index.md",
        "Lattices" => "lattice.md",
        "Finite lattices" => "finite_lattice.md",
        "Site ordering for MPS" => "mps_ordering.md",
        "Distances and neighbors" => "metric.md",
        "Operators and interactions" => "opsum.md",
        "Reading and writing files" => "io.md",
        "Plotting" => "plots.md",
        "Examples" => "examples.md",
    ],
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/awietek/Latlib.jl.git",
    devbranch = "main",
    push_preview = true,
)
