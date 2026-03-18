using TOML

@doc raw"""
    read_toml_interaction(filename::String; zero_based::Bool=true) -> OpSum

Reads a TOML file (such as the ones produced by `write_toml`) and returns the entries
listed under `Interactions` as an `OpSum` object.

Each interaction entry in the TOML file is expected to have the format:
```toml
Interactions = [
  ['coupling', 'type', site1, site2],
  ...
]
```
where 'coupling' and 'type' are strings, and `site1` and `site2` are integers.
The sites returned as OpSum are ALWAYS CONVERTED TO 1-BASED INDEXING!

# Arguments
- `filename::String`: Path to the TOML file to be read.

# Keyword arguments
- `zero_based::Bool=true`: If `true` (default), the site indices in the file
  are assumed to be 0-based and are converted to 1-based for the returned `OpSum`.
  Set to `false` if the file already uses 1-based indexing.
  Errors if index 0 is encountered while "zero_based" in set to false.
"""
function read_toml_interaction(filename::String; zero_based::Bool=true) :: OpSum
    data = TOML.parsefile(filename)
    if !haskey(data, "Interactions")
        error("TOML file does not contain an 'Interactions' section.")
    end

    # for 0-based file add 1, otherwise leave untouched
    offset = zero_based ? 1 : 0

    opsum = OpSum()
    for entry in data["Interactions"]
        cpl  = string(entry[1])
        type = string(entry[2])
        s1   = Int64(entry[3]) + offset
        s2   = Int64(entry[4]) + offset
        if !zero_based
            if s1 == 0 || s2 == 0
                error("0-based entries detected in `Interactions` section of TOML file $filename, which was requested to be read with keyword `zero_based=false`. It looks like the file does use 0-based indexing, so please read using `zero_based=true`.")
            end
        end
        opsum += Op(type, cpl, [s1, s2])
    end

    return opsum
end
