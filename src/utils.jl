import Base

# have a function that checks if a number is close to an integer within a certain tolerance
is_whole(x; atol=1e-8) = abs(x - Base.round(x)) ≤ atol

function get_latlib_version()
    # open Project.toml and read the version number
    open(joinpath(@__DIR__, "..", "Project.toml"), "r") do f
        for line in eachline(f)
            if startswith(line, "version")
                return split(line, "=")[2] |> strip
            end
        end
    end
end