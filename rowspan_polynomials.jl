using Catalyst, LinearAlgebra, InvertedIndices, SymPy, Symbolics

function rowspan_polynomials(reaction_network, complex_id_vec)
    M = complexstoichmat(reaction_network)*laplacianmat(reaction_network);

    n = size(M, 1);
    m = size(M, 2);

    N = M[:, Not(complex_id_vec)]

    dm = Int(size(M,1) - rank(SymPy.Sym.(symbolics_to_sympy.(M))))
    dn = Int(size(N,1) - rank(SymPy.Sym.(symbolics_to_sympy.(N))))

    if dn == dm 
        println("No rowspan polynomials exist involving the given complexes")
        return
    else

    e = Matrix{Int}(I, m, m)

    M_a = vcat(M, e[complex_id_vec, :])

    #Find the left nullspace of M_a 

    all_symbolics = unique(reduce(vcat, Symbolics.get_variables.(M_a)))

    converted_Ma = symbolics_to_sympy.(M_a) #convert to SymPy matrix

    Sp_Ma = SymPy.Sym.(converted_Ma)

    Sp_ns = reduce(hcat, transpose(Sp_Ma).nullspace())

    converted_ns = sympy_to_symbolics.(Sp_ns, Ref(all_symbolics))

    NlMa = Num.(converted_ns)

    # NlMa = nullspace(transpose(M_a))

    # setpoint = Symbolics.simplify(-NlMa[n+1,dm+1])
    num_binoms = length(complex_id_vec) - 1
    b = -NlMa[1:n,end-num_binoms+1:end]

    binoms = Symbolics.simplify.(transpose(b)*M)
    invariant = Symbolics.simplify.(sum(binoms, dims=1))
     
    if any(isequal.(invariant[complex_id_vec], 0))
        @warn "No linear invariant exists involving ALL the given complexes"
    end

    return b, invariant
    end
end