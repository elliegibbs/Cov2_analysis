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

# ╔═╡ fa4df7e0-86c5-488a-a9c2-fceaa93e5b81
histogram(keep_longer("../data/NCBI_sequences_download.fasta"),title = "Length of COVID Sequences (excluding <25k)", xlabel = "Sequence Length", ylabel = "Number of Sequences", leg = false)

# ╔═╡ 5e978c89-607b-49b9-bcc1-7b215eb40ca7
function delete_short_header
	too_short = length(sequences[2]) <= 25000
	findall(too_short)
end

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
# ╠═fa4df7e0-86c5-488a-a9c2-fceaa93e5b81
# ╠═5e978c89-607b-49b9-bcc1-7b215eb40ca7
