### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 0e972e5c-6b99-4833-9efe-7ecea0435dc7
begin
	using Pkg
	Pkg.develop(path="../../FinalBisc195.jl")
	using FinalBisc195
end

# ╔═╡ 5172e46e-885e-41ae-98f3-9190d1c7d193
using Plots

# ╔═╡ b0ce5162-e8d3-11eb-1b46-e71e1b2a5dda
md"# First Draft Covid Analysis Final"

# ╔═╡ 1a8310bb-8959-449f-b8af-a99709fcc0e9
md"### Analysis #1: Characterizing COVID Proteins"

# ╔═╡ d83d0535-da41-47bd-90f7-0394a9a2feea
# This is not something you should do, I'm doing it so I can run your other code
# You should fix your package code to handle ambiguous bases.
open("../data/fixseqs.fasta", "w") do io
	for line in eachline("../data/NCBI_sequences_download.fasta")
		if startswith(line, ">")
			println(io, line)
		else
			println(io, join([c in "ATGC" ? c : "N" for c in line]))
		end
	end
end

# ╔═╡ 0311e66c-8bae-459b-801c-4d5cbc6a9206
# TODO: previous version gets "invalid base R"
sequences = parse_fasta("../data/fixseqs.fasta") 

# ╔═╡ c04785aa-5b63-4d3d-846b-26c736d086b6
md"""
## To Do.

TODO: Use Markdown strings instead of comments

- Like this
- etc ...
"""

# ╔═╡ 08f99a4e-4d19-4795-9448-7e145d7e7271
#Here is a dictionary of codons translated into amino acids. This took time but it was mostly grunt work. 

# ╔═╡ 98708171-7c0d-4cc1-a346-5ff44fca4f32
function protein_dict(codon)
	proteins = Dict("GCT" => 'A',
					"GCC" => 'A',
					"GCA" => 'A',
					"GCG" => 'A',
					"GCT" => 'R',
					"GCG" => 'R',
					"CGA" => 'R',
					"CGG" => 'R',
					"AGA" => 'R',
					"AGG" => 'R',
					"CGT" => 'R',
					"CGC" => 'R',
					"AAT" => 'N',
					"AAC" => 'N',
					"GAT" => 'D',
					"GAC" => 'D',
					"TGT" => 'C',
					"TGC" => 'C',
					"CAA" => 'Q',
					"CAG" => 'Q',
					"GAA" => 'E',
					"GAG" => 'E',
					"GGT" => 'G',
					"GGC" => 'G',
					"GGA" => 'G',
					"GGG" => 'G',
					"CAT" => 'H',
					"CAC" => 'H',
					"ATT" => 'I',
					"ATC" => 'I',
					"ATA" => 'I',
					"ATG" => "*", # "*" means START... It also means "Methionine"
					"TTA" => 'L',
					"TTG" => 'L',
					"CTT" => 'L',
					"CTC" => 'L',
					"CTA" => 'L',
					"CTG" => 'L',
					"AAA" => 'K',
					"AAG" => 'K',
					"ATG" => 'M', # This overwrites the other one
					"TTT" => 'F',
					"TTC" => 'F',
					"CCT" => 'P',
					"CCC" => 'P',
					"CCA" => 'P',
					"CCG" => 'P',
					"TCT" => 'S',
					"TCC" => 'S',
					"TCA" => 'S',
					"TCG" => 'S',
					"AGT" => 'S',
					"AGC" => 'S',
					"ACT" => 'T',
					"ACC" => 'T',
					"ACA" => 'T',
 					"ACG" => 'T',
					"TGG" => 'W',
					"TAT" => 'Y',
					"TAC" => 'Y',
					"GTT" => 'V',
					"GTC" => 'V',
					"GTA" => 'V',
					"GTG" => 'V',
					"TAA" => "**", #** means STOP
					"TGA" => "**",
					"TAG" => "**")
		!(codon in keys(proteins)) && error("Invalid codon $codon")
		return proteins[codon]
end

# ╔═╡ e3644d03-da45-48ec-a9b2-b1212add7d5e
#This function takes a sequence and divides it into groups of three (codons) and then  runs each codon through the dictionary. It was very difficult to write, especailly the codon = SubString..etc line. I had to use Google to get the formula necessary to parse the String. After the codons have been parsed and translated, they are put into another array then joined into a string of amino acids. 

# ╔═╡ 3d4ab822-a2ec-4851-b449-f92b42e7bb06
md"""
Re: comment ☝

I'm not super clear why you did it this way instead of the way that we pulled out kmers before, eg `sequence[i:i+2]`. 
You can also iterate through ranges by threes, eg `1:3:stopping_point`
"""

# ╔═╡ 57c894ac-6e63-4f35-bc33-415dba0ca36d
some_seq = "ATGGGGCTGAGAGA"

# ╔═╡ 8089a1b2-3f26-497d-96ac-57daf2f290c6
begin
	codons = String[]
	for i in 1:3:length(some_seq)-2
		push!(codons, some_seq[i:i+2])
	end
	codons
end

# ╔═╡ da3ceb6b-9640-4df6-af58-57599679eb9a
# Do you know that your sequences always start at position 1?
# If you want to do this, you need to add logic to find the start codon
# Also note that your genomes actually contain multiple protein coding genes, not just one
# you could also just download the protein sequences from NCBI directly
function convert_protein(sequence)
	protein_array = []
	codon_array =[]
	sequence_length = length(sequence)
	third_sequence_length = trunc(Int,((sequence_length)/3))
	for i in 1:third_sequence_length
		if i <= (length(sequence)-2)
			codon = SubString(sequence, 3(i-1)+1,3(i-1)+3)
			push!(codon_array, codon)
			i = i+3
		end
	end
	
	for i in 1:length(codon_array)
		codon = protein_dict(codon_array[i])
		push!(protein_array, codon)
		i=i+1
	end
	protein_string = join(protein_array)
	return protein_string
end 
	


# ╔═╡ e67ecd4b-a7d3-4adc-ad7d-27f9a23279fe
md"""
Looks like this sequence stops after ~10 amino acids?

Pretty sure this translation is out of frame
"""

# ╔═╡ c1267b6d-6012-4e49-b3dd-e3414d4da611
protein_one= convert_protein(sequences[2][1])

# ╔═╡ dbd2e3a8-4500-4a78-b82c-f4295920d9d6
#This function will be used for analysis of the amino acids. It classifies them as basic, acidic, polar, and nonpoalr. However, it can only do this for one protein sequence at a time. 

# ╔═╡ 0533ea5a-10c6-4fb0-821b-339de0c5a3b1
function classify_proteins(protein_sequence)
	basic = "KHR"
	acidic = "DE"
	polar = "GSTCYNQ"
	nonpolar = "AVLIPFWM"

	
	basic_amino = 0
	acidic_amino = 0
	polar_amino = 0
	nonpolar_amino = 0
	
	for aminoacid in protein_sequence
		if occursin(aminoacid, basic)
			basic_amino= basic_amino+1
		elseif occursin(aminoacid, acidic)
			acidic_amino = acidic_amino +1
		elseif occursin(aminoacid, polar)
			polar_amino = polar_amino+1
		elseif occursin(aminoacid, nonpolar)
			nonpolar_amino = nonpolar_amino+1
		end
	end
	return basic_amino, acidic_amino, polar_amino, nonpolar_amino
end

# ╔═╡ 4aaddc38-c25e-465b-ae20-3e36bdfc3ea7
classify_proteins(protein_one)

# ╔═╡ 4f7a893a-1be7-45b2-a2f9-703f93eaf513
#Now, I need to make a function that takes each analysis and puts it into one array so that I can make a histogram that will compare them. 

# ╔═╡ 6002935f-859b-4f81-ba1a-36736c19537b
function compare_proteins(path)
	header, seqs = parse_fasta(path)
	analyzed = []
	records = count_records(path)
	for i in seqs
			protein = convert_protein(seqs[i])
			@info protein
			analyze_protein = classify_proteins(protein)
			@info analyze_protein
			push!(analyzed, analyze_protein)
		end
	return analyzed
end

# ╔═╡ 5bcf2bc9-71c7-48bb-8642-451a6e5d7b79
compare_proteins("data/cov2_genomes.fasta") # what's this file? Not mentioned in `data.md`

# ╔═╡ 7372f03c-8d14-4f52-8d61-9bc0719e4622
# where is `seq_lengths` defined? It should be in your package
histogram(seq_lengths("data/NCBI_sequences_download.fasta"),title = "Length of COVID Sequences", xlabel = "Sequence Length", ylabel = "Number of Sequences", leg = false)

# ╔═╡ Cell order:
# ╟─b0ce5162-e8d3-11eb-1b46-e71e1b2a5dda
# ╟─1a8310bb-8959-449f-b8af-a99709fcc0e9
# ╠═0e972e5c-6b99-4833-9efe-7ecea0435dc7
# ╠═d83d0535-da41-47bd-90f7-0394a9a2feea
# ╠═0311e66c-8bae-459b-801c-4d5cbc6a9206
# ╠═c04785aa-5b63-4d3d-846b-26c736d086b6
# ╠═08f99a4e-4d19-4795-9448-7e145d7e7271
# ╠═98708171-7c0d-4cc1-a346-5ff44fca4f32
# ╠═e3644d03-da45-48ec-a9b2-b1212add7d5e
# ╠═3d4ab822-a2ec-4851-b449-f92b42e7bb06
# ╟─57c894ac-6e63-4f35-bc33-415dba0ca36d
# ╠═8089a1b2-3f26-497d-96ac-57daf2f290c6
# ╠═da3ceb6b-9640-4df6-af58-57599679eb9a
# ╟─e67ecd4b-a7d3-4adc-ad7d-27f9a23279fe
# ╠═c1267b6d-6012-4e49-b3dd-e3414d4da611
# ╠═dbd2e3a8-4500-4a78-b82c-f4295920d9d6
# ╠═0533ea5a-10c6-4fb0-821b-339de0c5a3b1
# ╠═4aaddc38-c25e-465b-ae20-3e36bdfc3ea7
# ╠═4f7a893a-1be7-45b2-a2f9-703f93eaf513
# ╠═6002935f-859b-4f81-ba1a-36736c19537b
# ╠═5bcf2bc9-71c7-48bb-8642-451a6e5d7b79
# ╠═5172e46e-885e-41ae-98f3-9190d1c7d193
# ╠═7372f03c-8d14-4f52-8d61-9bc0719e4622
