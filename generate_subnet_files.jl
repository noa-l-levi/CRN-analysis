function generate_subnet_files(net_name::Symbol)
    str1 = string(net_name)
    include(str1 * ".jl")

    # Now get the actual ReactionSystem object by name
    sys = eval(net_name) 

    subnet_indices = independent_decomposition(sys)

    for (i, idxs) in enumerate(subnet_indices)
        println("Subnetwork $i: ", idxs)
        filename = "$(str1)_SN$(i).jl"

        # Build code strings
        include_line = "include(\"$str1.jl\")"
        using_lines = "using Catalyst\nusing ModelingToolkit"
        rs_line = "@named subnet = ReactionSystem(reactions($str1)[$(repr(idxs))], default_t())"
        complete_line = "$(str1)_SN$(i) = complete(subnet)"

        file_content = """
        $include_line
        $using_lines

        $rs_line
        $complete_line
        """

        open(filename, "w") do f
            write(f, file_content)
        end
    end
end