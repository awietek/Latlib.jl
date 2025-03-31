using Documenter, Latlib
#using LiveServer

makedocs(sitename="Latlib Documentation")

#LiveServer.serve(dir="/Users/soares/.julia/dev/Latlib/docs/build", port=8000)

deploydocs(
    repo = "github.com/awietek/Latlib.jl.git",
)
