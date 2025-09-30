using Catalyst, Oscar, GAP, RowEchelon

function get_Oscar_SteadyStatesystem(rn)
    #Set up the equations
    #INPUT: rn
    #OUTPUT: system of polynomial equations
        
        N = netstoichmat(rn) #stoichiometric matrix
        Y = substoichmat(rn) #reactant matrix, (columns correspond to exponent vectors of the monomials)
        W = Oscar.rref(matrix(QQ,conservationlaws(N)))[2]
        n,r = size(N) #number of species and reactions
        d = size(W)[1] #number of conservation laws
        
        
        #Set up the polynomial system
        global T, k,c = polynomial_ring(QQ, "k"=>1:r,"c"=>1:d;internal_ordering=:deglex)
        global K = fraction_field(T)
        global B, x = polynomial_ring(K, "x"=>1:n;internal_ordering=:deglex)
        # polynomialSystem = Array{fmpq_mpoly, 1}()
        polynomialSystem = Vector{elem_type(B)}()
         
        #array for the equation system
        F = []
        
        #Create the monomial with exponent Y[:,j]
        monomials = []
    
        for j in 1:r
            mono = B(1)
            for i in 1:n
                mono = mono*x[i]^Y[i,j]
            end
            push!(monomials,mono)
        end
    
        #Create the steady state equations
        for i in 1:n
            f = B(0)
            for j in 1:r
                #Create the monomial with exponent Y[:,j]
                f = f + N[i,j]*k[j]*monomials[j]
            end
            push!(F,f)
        end
        
      
        return F
    end