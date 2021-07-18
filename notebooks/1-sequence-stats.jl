### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ d77a3d79-f377-4540-8b63-0af54f961db9
md"# Assignment 8"

# ╔═╡ fda59124-5dc1-488e-b7fd-b01b854d8dba
# Calculate the mean sequence length and standard deviation of the CoV2 geneomes

# ╔═╡ 51a063ff-6d72-453c-91b7-8882235ff8a6
function mean_cov2(seqs)
    lengths = []
    for seq in seqs
        push!(lengths, length(seq))
    end
    mean_cov = sum(lengths) / length(headers)
    return mean_cov
end

# ╔═╡ 0a673e6a-79b5-424c-858c-55f446cfb15c
function stand_dev(seqs)
    lengths = []
    for seq in seqs
        push!(lengths, length(seq))
    end
    samples_mean = mean_cov2(seqs)
    samples_size = length(headers)
    samples = map(x -> (x - samples_mean)^2, lengths)
    samples_sum = sum(lengths)
    samples_std = sqrt(samples_sum / (samples_size - 1))
    return samples_std
end

# ╔═╡ 8bb7c8a1-fe79-497b-8d13-41cd968dcb84
#Calculate the gc content of a sequence

# ╔═╡ cc3c8cc1-6e2f-462a-8b9d-9fc31eb8b3c3
function gc_content(sequence)
    upsequence = uppercase(sequence)
        seqlength = length(upsequence)
        gs = count(==('G'), upsequence)
        cs = count(==('C'), upsequence)
        return (gs + cs) / seqlength 
end

# ╔═╡ 6ca306e0-15f9-4c76-99de-c25cd7c64000
# Calculate the mean and standard deviation of GC content of the CoV2 genomes

# ╔═╡ cb86bbf9-7102-4724-b300-f6487bd9eed4
function mean_gc(path)
    headers, bodies = parse_fasta(path)
    gc_array=[]
    for body in bodies
        each_gc = gc_content(body)
        push!(gc_array, each_gc)
    end
    array_sum = sum(gc_array)
    gc_mean = array_sum / count_records(path)
    return gc_mean
end

# ╔═╡ b61086cd-f163-4a04-ae82-b1d2c3bd3212
function gc_stand_dev(path)
    headers, bodies = parse_fasta(path)
    gc_array=[]
    for body in bodies
        each_gc = gc_content(body)
        push!(gc_array, each_gc)
    end
    samples_mean = mean_gc(path)
    samples_size = count_records(path)
    samples = map(x -> (x - samples_mean)^2, gc_array)
    samples_sum = sum(samples)
    samples_std = sqrt(samples_sum / (samples_size - 1))
    return samples_std
end

# ╔═╡ Cell order:
# ╟─d77a3d79-f377-4540-8b63-0af54f961db9
# ╠═fda59124-5dc1-488e-b7fd-b01b854d8dba
# ╠═51a063ff-6d72-453c-91b7-8882235ff8a6
# ╠═0a673e6a-79b5-424c-858c-55f446cfb15c
# ╠═8bb7c8a1-fe79-497b-8d13-41cd968dcb84
# ╠═cc3c8cc1-6e2f-462a-8b9d-9fc31eb8b3c3
# ╠═6ca306e0-15f9-4c76-99de-c25cd7c64000
# ╠═cb86bbf9-7102-4724-b300-f6487bd9eed4
# ╠═b61086cd-f163-4a04-ae82-b1d2c3bd3212
