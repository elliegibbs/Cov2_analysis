### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ ad0defb1-f8e1-4c5c-ac80-337199ceca47
using FinalBisc195

# ╔═╡ d45f2f9d-2c1a-4f4c-833d-c9c947fb065c
using Plots

# ╔═╡ d77a3d79-f377-4540-8b63-0af54f961db9
md"# Cov2 Analysis"

# ╔═╡ a0d3f8ae-4417-4d9c-81ff-77e6f66fb0c9
sequences = parse_fasta("../data/NCBI_sequences_download.fasta")

# ╔═╡ fda59124-5dc1-488e-b7fd-b01b854d8dba
# Calculate the mean sequence length and standard deviation of the CoV2 geneomes

# ╔═╡ 51d7bf33-baf7-432f-a071-07ffde9d92d0
mean_cov2(sequences[2])

# ╔═╡ dd8c5549-03f1-42b4-8b3f-552b83f61526
stand_dev(sequences[2])

# ╔═╡ 6ca306e0-15f9-4c76-99de-c25cd7c64000
# Calculate the mean and standard deviation of GC content of the CoV2 genomes

# ╔═╡ 29c771a8-5436-4e80-bd98-f29ae57555e9
mean_gc("../data/NCBI_sequences_download.fasta")

# ╔═╡ d7be8ded-9ed2-4913-a181-de6066666a86
gc_stand_dev("../data/NCBI_sequences_download.fasta")

# ╔═╡ ebed478c-7108-4523-9ee6-4f6d6dac6ebb
#find the length of the shortest and longest sequences (shortest, longest)

# ╔═╡ 804fa081-9346-4408-b695-f8473959e971
shortest_and_longest("../data/NCBI_sequences_download.fasta")

# ╔═╡ 908c6eea-d397-4395-b7ac-649499a0840b
histogram(seq_lengths("../data/NCBI_sequences_download.fasta"),title = "Length of COVID Sequences", xlabel = "Sequence Length", ylabel = "Number of Sequences", leg = false)

# ╔═╡ a6ac558c-c4d9-4af9-89a1-76dac03f33c5
function keep_longer(path)
    headers, bodies = parse_fasta(path)
    length_array = seq_lengths(path)
	good_bodies_lengths = []
	good_heads = []
	
	
    good_bodies_index = findall(x ->x >= 25000, length_array)
	
	

	for i in good_bodies_index
		push!(good_bodies_lengths, length_array[i])
		push!(good_heads, headers[good_bodies_index])
	end
	
	
	return good_bodies_lengths
	
	
    
end

# ╔═╡ e448be6e-a547-48de-ad25-318a06a2a662
keep_longer("../data/NCBI_sequences_download.fasta")

# ╔═╡ fa4df7e0-86c5-488a-a9c2-fceaa93e5b81
histogram(keep_longer("../data/NCBI_sequences_download.fasta"),title = "Length of COVID Sequences (excluding <25k)", xlabel = "Sequence Length", ylabel = "Number of Sequences", leg = false)

# ╔═╡ 11b509db-48b2-4d3a-bb01-b1007d230c3e
unique_kmer_count(sequences[2][1], 3)

# ╔═╡ b77034dd-b4b2-4cc5-b0c6-03a0d0a2718e
function  kmercount(sequence, k)
    1 <= k <= length(sequence) || error("k must be a positive integer less than the length of the sequence")
    kmers = Dict() # initialize dictionary

    stopindex = length(sequence)-(k-1)

    for i in 1:stopindex
        kmer = sequence[i:i+k-1] # Change to index the sequence from i to i+k-1
        kmer = normalizeDNA(kmer) 
        indict = haskey(kmers, kmer)
        if indict == true
            println(true)
            kmers[kmer]=kmers[kmer]+1
        else 
            kmers[kmer]=1
            println(false)
        end
    end

        #   if this kmer is a key the dictionary
        #       add 1 to the value referenced by that kmer
        #   otherwise
        #       make a new entry in the dictionary with a value of 1

    return kmers
end

# ╔═╡ 6ec2027d-dcd3-4c9e-92c3-4a3694e3fb1c
function kmer_distance(kmer1, kmer2)
	length_of_intersection= length(intersect(kmer1,kmer2))
	length_of_union=length(union(kmer1,kmer2))
	distance_metric = 1-(length_of_intersection/length_of_union)
end

# ╔═╡ 26a48a4a-2e59-4247-859b-c45858588fa2
first_kmer=kmercount(sequences[2][1],7)



# ╔═╡ 8287bb96-6fbf-4abf-9241-bac28e4af855
second_kmer=kmercount(sequences[2][2],7)

# ╔═╡ 08488f76-6b35-4875-bffe-25130925a796
kmer_distance(first_kmer, second_kmer)

# ╔═╡ Cell order:
# ╟─d77a3d79-f377-4540-8b63-0af54f961db9
# ╠═ad0defb1-f8e1-4c5c-ac80-337199ceca47
# ╠═a0d3f8ae-4417-4d9c-81ff-77e6f66fb0c9
# ╠═fda59124-5dc1-488e-b7fd-b01b854d8dba
# ╠═51d7bf33-baf7-432f-a071-07ffde9d92d0
# ╠═dd8c5549-03f1-42b4-8b3f-552b83f61526
# ╠═6ca306e0-15f9-4c76-99de-c25cd7c64000
# ╠═29c771a8-5436-4e80-bd98-f29ae57555e9
# ╠═d7be8ded-9ed2-4913-a181-de6066666a86
# ╠═ebed478c-7108-4523-9ee6-4f6d6dac6ebb
# ╠═804fa081-9346-4408-b695-f8473959e971
# ╠═d45f2f9d-2c1a-4f4c-833d-c9c947fb065c
# ╠═908c6eea-d397-4395-b7ac-649499a0840b
# ╠═a6ac558c-c4d9-4af9-89a1-76dac03f33c5
# ╠═e448be6e-a547-48de-ad25-318a06a2a662
# ╠═fa4df7e0-86c5-488a-a9c2-fceaa93e5b81
# ╠═11b509db-48b2-4d3a-bb01-b1007d230c3e
# ╠═b77034dd-b4b2-4cc5-b0c6-03a0d0a2718e
# ╠═6ec2027d-dcd3-4c9e-92c3-4a3694e3fb1c
# ╠═26a48a4a-2e59-4247-859b-c45858588fa2
# ╠═8287bb96-6fbf-4abf-9241-bac28e4af855
# ╠═08488f76-6b35-4875-bffe-25130925a796
