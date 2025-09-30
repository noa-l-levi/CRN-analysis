# loading the required packages
using RowEchelon, Catalyst, Graphs, Combinatorics, MetaGraphs, ModelingToolkit

# independent decomposition function
function independent_decomposition(reaction_network)

    # compute the deficiency
    reactioncomplexes(reaction_network);
    linkageclasses(reaction_network);
    numspecies(reaction_network);
    N = netstoichmat(reaction_network);
    δ = deficiency(reaction_network);

    # get the transpose of the stoichiometric matrix
    R = transpose(N);

    # compute a maximal linearly independent set of vectors in R
    pivots = rref_with_pivots(transpose(R));

    basis = R[pivots[2],:];

    # initialise an undirected graph G
    G = SimpleGraph();
    G = MetaGraph(G);

    # add vertices to G which are the vectors in basis
    for i = 1:length(pivots[2]);
        Graphs.add_vertex!(G);
        st = string("R",pivots[2][i]);
        set_prop!(G, i, :name, st);
    end;

    # initialise matrix of linear combinations (#rxns x #basis elements)
    linear_combo = zeros(size(N)[2], length(pivots[2]));

    # get matrix which shows the linear combinations of basis rxns which 
    # give the nonbasis rxns (basis rxns will have a row of zeros)
    for i = 1:size(N)[2];
        if !in.(i, Ref(pivots[2]));

            linear_combo[i,:] = transpose(basis)\R[i,:];
        end;
    end;

    linear_combo = round.(linear_combo, digits=1); 

    # get the reactions that are linear combinations of at least 2 basis reactions
    # these give us the edges
    abs_linear_combo = abs.(linear_combo);
    get_edges = first.(Tuple.(findall(x -> x > 1, sum(abs_linear_combo, dims=2)))); 

    #  initialize an array for sets of vertices that will form the edges
    vertex_set = Array{}[];
        
    # identify which vertices form edges in each reaction: 
    # get those with non-zero coefficients in the linear combinations
    for i = 1:length(get_edges);
        push!(vertex_set, [first.(Tuple.(findall(!iszero, linear_combo[get_edges[i], :])))]);
    end;
        
    # initialize the edge set
    edges = Vector{}[];
        
    # get all possible combinations (not permutations) of the reactions 
    # involved in the linear combinations
    for i = 1:length(vertex_set);
        append!(edges, [collect(combinations(vertex_set[i][1], 2))]);
    end;
        
    #  get just the unique edges (if edges is non-empty)
    if isempty(edges) != true;
        edges = unique(reduce(vcat,edges));
    end;
        
    # add these edges to graph G
    for i = 1:length(edges);
        Graphs.add_edge!(G, edges[i][1], edges[i][2]); 
    end;

    # determine whether a finer decomposition exists, if not, end here
    if Graphs.is_connected(G) == true
        println("No finer decomposition exists")
        return
    end
    
    if Graphs.is_connected(G) != true
        println("A finer network decomposition exists")
    end

    # determine which vertices are connected
    components = connected_components(G);
        
    # determine the number of connected components of G 
    # this is the number of partitions R will be decomposed to
    println("The number of independent subnetworks is:")
    num_components = size(components,1)
    display(num_components)

    # initialize the vector of partitions
    P = [];
        
    # basis vectors: assign them first into their respective partition 
    for i = 1:size(components)[1];
        push!(P,pivots[2][components[i][1:size(components[i])[1]]]);
    end;
    
    # nonbasis vectors: they go to the same partition as the basis vectors 
    # that form their linear combination
    for i = 1:size(P,1);
        for j = 1:size(P[i],1);
            
            # get the column number representing the basis vectors in 'linear_combo'
            col = findall(x -> x == P[i][j], pivots[2])[1];
            
            # check which reactions used a particular basis vector and assign them 
            # to their respective partition
            append!(P[i], findall(x -> x != 0, linear_combo[:, col])); 
        end;
    end;
    
    # get only unique elements in each partition
    for i = 1:size(P,1);
        P[i] = unique(P[i]);
    end;
    return(P)


    # check that all reactions have been assigned a partition
    if length(reduce(vcat,P)) == size(R, 1) 
        println("No reactions have been omitted in the partition") 
        # return
    end
    
    # if reactions are missing, repeat from the linear combination step
    if length(reduce(vcat,P)) != size(R, 1) 
        println("Warning: some reactions have been omitted in the partition") 
        missing = size(R, 1) - length(reduce(vcat,P))
        display(missing)
   
    end

    
end
