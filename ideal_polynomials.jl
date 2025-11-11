using Symbolics, Catalyst
using Catalyst: unknowns 
using Oscar
using Oscar: groebner_basis
using InvertedIndices: Not

function ideal_polynomials(reaction_network, select_vars) 

    function polynomial_ideal(eqn_rhs, vars, params)
        # Convert Catalyst symbolic variables into Julia Symbol types
        varnames = tosymbol.(vars, escape=false)
        paramnames = tosymbol.(params)
        
        # Create polynomial ring in Oscar
        CC, oscar_coeffs = polynomial_ring(QQ, paramnames)
        ff = fraction_field(CC)
        RR, oscar_vars = polynomial_ring(ff, varnames)
        
        # Map Catalyst variables to Oscar variables
        cat_var_params = [vars; params]
        oscar_var_params = [oscar_vars; oscar_coeffs]
        cat_to_oscar = Dict(cat => oscar for (cat, oscar) in zip(cat_var_params, oscar_var_params))
        oscar_to_cat = Dict((oscar => cat) for (cat, oscar) in cat_to_oscar)
        
        # Build Oscar polynomial by substituting oscar vars in Catalyst equations RHS (Right Hand Sides)
        polys = map(eqn_rhs) do rhs
            if rhs isa Number  # is a constant e.g. zero
                RR(rhs)
            else
                Symbolics.substitute(rhs, cat_to_oscar)
            end
        end
        ideal(RR, polys)
    end

    function polynomial_ideal(odesys::ODESystem)
        vars = unknowns(odesys)
        params = parameters(odesys)
        eqn_rhs = [eq.rhs for eq in equations(odesys)]
        polynomial_ideal(eqn_rhs, vars, params)
    end


    function to_symbolic_polynomial(poly::MPolyRingElem)
        # create new symbolic vars from ring of poly
        var_syms = Symbolics.variable.(gens(parent(poly)))
        coeff_syms = Symbolics.variable.(gens(coefficient_ring(poly)))
    
        # coefficients
        coeff_terms = Oscar.coefficients(poly)
        coeffs = zeros(Num, length(coeff_terms))
        for (i, term) in enumerate(coeff_terms)
            cf_and_es = coefficients_and_exponents(term.num)
            coeff_polyn = sum(Int(cf) * prod(coeff_syms .^ es) for (cf, es) in cf_and_es)
            coeffs[i] = coeff_polyn
        end

        # polynomial variables
        exp_vecs = collect(Oscar.exponents(poly))
        xs = [prod(var_syms .^ e_vec) for e_vec in exp_vecs]
        coeffs' * xs
    end

    # function lift_symbolics(f, I, xs)
    #     A = coordinates(f, I)
    #     B = zeros(Num, length(A))
    #     for (i, aa) in enumerate(A)
    #         if length(aa.coeffs) == 0
    #             continue
    #         else
    #             # Process all coefficients in the polynomial
    #             for ff_coeff in aa.coeffs
    #                 cf_and_es = coefficients_and_exponents(ff_coeff.num)
    #                 polyn = sum(Int(cf) * prod(xs .^ es) for (cf, es) in cf_and_es)
    #                 B[i] += polyn  # Accumulate all coefficient contributions
    #             end
    #         end
    #     end
    #     return B
    # end

    function lift_symbolics(f, I, xs)
        A = coordinates(f, I)
        B = zeros(Num, length(A))
        for (i, aa) in enumerate(A)
            if length(aa.coeffs) == 0
                continue
            else
                @show aa.coeffs
                ff_coeff = aa.coeffs[1]
                cf_and_es = coefficients_and_exponents(ff_coeff.num)
                polyn = sum(Int(cf) * prod(xs .^ es) for (cf, es) in cf_and_es)
                B[i] = polyn
            end
        end
        return B
    end


    odesys = convert(ODESystem, reaction_network, combinatoric_ratelaws=false)
    I = polynomial_ideal(odesys)
    xs = gens(base_ring(I))

    # Confirm selected species
    println("Selected species: ", string.(xs[select_vars]))
    
    elimination_ordering = degrevlex(xs[Not(select_vars)]) * degrevlex(xs[select_vars]) 

    # Confirm elimination ordering
    if is_elimination_ordering(elimination_ordering, xs[Not(select_vars)]) == false
        error("The specified ordering is not an elimination ordering for the selected variables.")
    end
   
    G = Oscar.groebner_basis(I; ordering=elimination_ordering)

    elimination_ideal =[]
    for k in 1:length(G)
        var_indices = [i for i in 1:length(xs) if Oscar.degree(G[k], i) > 0]
        unselected_indices = filter(x -> !(x in select_vars), var_indices)
        if isempty(unselected_indices)
            push!(elimination_ideal, G[k])
            println("Elimination ideal element ", k, ":", G[k])
        end
    end

    lift_output = []
    elim_ideal_polys = []
    for k in 1:length(elimination_ideal)
        params = parameters(reaction_network)
        # sym_coeffs = lift_symbolics(G[k], I, params)
        sym_coeffs = lift_symbolics(elimination_ideal[k], I, params)

        eqn_lhs = [eqn.lhs for eqn in equations(odesys)]
        eqn_rhs = [eqn.rhs for eqn in equations(odesys)]

        poly_comb_k = dot(eqn_lhs, sym_coeffs) |> Symbolics.expand 
        push!(lift_output, poly_comb_k)

        poly_k = dot(eqn_rhs, sym_coeffs) |> Symbolics.expand
        push!(elim_ideal_polys, poly_k)
    end

    return G, elim_ideal_polys, lift_output

end