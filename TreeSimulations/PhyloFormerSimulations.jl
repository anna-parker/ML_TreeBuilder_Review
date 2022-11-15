
using ARGTools
using TreeTools
using Combinatorics
using DataFrames
using TreeAlgs.Evolve
using ArgParse
using Random
using CSV

function get_r(ρ, n, N, simtype::Symbol)
    if simtype == :kingman
        return ρ * n / N
    elseif simtype == :yule
        return ρ / N
    elseif simtype == :flu
    	return ρ * n^0.2 / N
    else
        @error "Unrecognized `simtype`."
    end
end

"""
get_trees(no_trees, no_lineages; remove=false, debug=false, c=0.75, ρ = 0.05
returns trees as well as true MCCs
simulate a total of `no_trees` trees using the ARGTools package, 
specifying the `lineage_no` determines the number of nodes in each tree

ρ - Reassortment rate scaled to coalescence rate
s - Fraction of leaves that are not sampled at time 0
"""
function get_trees(no_trees, no_lineages; ρ = 0.05, N = 10_000, simtype = :flu, s=0.0)
    # Parameters of the ARG simulation
    N = N # pop size
    nmax = no_lineages  # total number of lineages
    if s!=0.0 #if additional parameters should be added
        #n0 = max(1, round((1-s)*no_lineages))
        n0 = 0.25*nmax
        s = s*(nmax-1)*nmax^0.2 /(2*N) # add leaves before time 0 at a rate relative to coalescence
    else
        n0 = nmax
    end
    r = get_r(ρ, 0.75*nmax, N, simtype) # Absolute reassortment rate

    # Simulating the ARG
    arg = ARGTools.SimulateARG.simulate(N, r, n0; s, nmax, K=no_trees, simtype);
    # The trees for the 2 segments
    trees = ARGTools.trees_from_ARG(arg; node_data = TreeTools.MiscData);
    return trees, arg
end


function get_c(res_rate, rec; n=100, simtype=:flu)
    if simtype == :flu
        c_values = CSV.read(Base.pwd()*"/Performance_Influenza/influenza_c_values.txt", DataFrame)
        keys = c_values.rec_rate 
        if res_rate == 0.3
            values = c_values."res_0.3" 
        elseif res_rate == 0.35
            values = c_values."res_0.35" 
        elseif res_rate == 0.4
            values = c_values."res_0.4" 
        else
            throw(Error("not implemented"))
        end
        dict_ = Dict(zip(keys, values))
        return dict_[rec]
    elseif simtype == :kingman
        return 0.75
    end
end


function write_real_tree!(tree::Tree, μ::Number, tree_name::String, outfolder::AbstractString)
    #normalize tree by mutation rate
    for node in POT(tree)
        node.tau= node.tau * μ
    end

    write_newick(outfolder *"/"*tree_name *".nwk", tree)
end

"""
    simulate(n::Int, r::Number, outfolder::AbstractString; flu = true, debug=true)
Simulate a pair of ARG trees using sample size (`n`), and recombination rate ratio to coalescence (`r`),
mutate sequences along simulated trees (if `flu` use sequence lengths of HA and NA segments, else 1000).
"""
function simulate(n::Int, rec_rate::Number, outfolder::AbstractString; res_rate=0.3, simtype = :flu, s=0.0)
    ##set sequence length length and mutation rate
    r = 10^rec_rate
    L = 1000
    N = 10_000
    c = get_c(res_rate, rec_rate; n=n, simtype)
    μ = 4/(N*c*L*3)

    ##simulate trees and evolve segments
    trees, arg = get_trees(8, n; ρ = r, simtype, s)
    if !isdir(outfolder)
        mkdir(outfolder)
    end
    for i in 1:8
        trees[i].label = string(i)
        evolve!(trees[i], L, μ);
        write_seq2fasta(trees[i], string(i), outfolder, only_terminals=true, remove_0_mutations=false, write_date=true, year_rate=4.73e-3/μ);
        write_real_tree!(trees[i], μ, "true_"*string(i), outfolder)
    end

    delete_list = Vector{String}[]
    for tree in trees
        l = []
        for node in internals(tree)
            if !node.isroot
                if !node.data.dat["evolved"]
                    push!(l, node.label)
                end
            end
        end
        append!(delete_list, [l])
    end
    print("Average res rate: "*string(((n-1-sum([length(l) for l in delete_list])/length(delete_list))/(n-1)))*"\n")
    ## remove internal nodes (except the root) with no mutations between them and their children

    for i in 1:8
        to_be_deleted = delete_list[i]
        for node_label in to_be_deleted
            delete_node!(trees[i], node_label, ptau=true)
        end
    end

    ## write modified tree to nwk files
    mkpath(outfolder)
    for i in 1:8
        write_newick(outfolder *"/"*string(i)*".nwk", trees[i])
    end

end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--n"
            help = "sample size at time 0 (Int)"
            arg_type = Int
            default = 100
        "--rec"
            help = "recombination rate ratio to coalescence (Float)"
            arg_type = Float64
            default = 0.0
        "--o"
            help = "output directory"
            arg_type = String
            default = "results"
        "--simtype"
            help = "which simtype to use"
            arg_type = String
            default = "flu"
        "--res"
            help = "resolution rate"
            arg_type = Float64
            default = 0.3
        "--s"
            help= "rate to add nodes backwards in time"
            arg_type = Float64
            default = 0.0
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    n = parsed_args["n"]
    r = parsed_args["rec"]
    o = parsed_args["o"]
    simt = parsed_args["simtype"]
    if simt=="flu"
        simtype = :flu
    else
        simtype = :kingman
    end
    res = parsed_args["res"]
    s = parsed_args["s"]
    println("Simulating trees and sequences of sample size $n and recombination rate $r")
    simulate(n, r, o; res_rate=res, simtype, s)
end

main()