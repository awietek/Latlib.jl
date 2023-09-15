using Printf

color_cycle = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple",
               "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan"]

function plot(lat::Lattice, bonds::Vector{Bond}; annotate::Bool=true, zero_indexed::Bool=false)
    fig, ax = subplots()
    coords = coordinates(lat)

    if dimension(lat, true) == 2
        for (idx, c) in enumerate(eachcol(coords))
            ax.plot(c[1], c[2], "o", color="tab:gray", ms=10, zorder=10)
            if annotate
                if zero_indexed
                    ax.annotate(string(idx-1), (c[1], c[2]), zorder=20)
                else
                    ax.annotate(string(idx), (c[1], c[2]), zorder=20)
                end
            end
        end

        types_couplings = unique([(b.type, b.coupling) for b in bonds])
        for (tcidx, (t, c)) in enumerate(types_couplings)
            tcbonds = [b for b in bonds if t == b.type && c == b.coupling]
            label = @sprintf("%s %s", t, c)
            color = color_cycle[mod1(tcidx, length(color_cycle))]
            for b in tcbonds
                c1 = coords[:,b.sites[1]]
                c2 = coords[:,b.sites[2]]
                plt.plot([c1[1], c2[1]], [c1[2], c2[2]], color=color, label=label)
                label = nothing            
            end
        end
        
        ax.set_aspect("equal")
        ax.legend()
        show()
    end
end
