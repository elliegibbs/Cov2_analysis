### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 9cda8166-e37f-11eb-2e1b-bf74d988c460
md"### Final Assignment Analysis Plan"

# ╔═╡ e7ad56ce-336e-4c1b-9aae-c3a989f7174b
md"Analysis #1: Compare NW vs SW alignment"

# ╔═╡ a8a8631b-d883-480a-83e6-3015bffb2e3b
function NW(seq1, seq2)
	# use NW alignment function from lab
			function nwscore(base1::Char, base2::Char; match = 1, mismatch = -1, gap = 				-1)
			if base1 == base2
				return match
			else base1 != base2
				return mismatch
			end
		end

		function nwscore(base::Char, ::Nothing; match = 1, mismatch = -1, gap = -1)
			return gap
		end

		function nwscore(::Nothing, base::Char; match = 1, mismatch = -1, gap = -1)
			nwscore(base, nothing; gap)
		end

		function nwscore(::Nothing, ::Nothing)
			throw(ArgumentError("Score for two gaps is not defined"))
		end
			function nwsetupmatrix(s1, s2; gap=-1)
			setupmatrix = zeros(Int, length(s1)+1, length(s2)+1)
			for j in 1:length(s2)+1
				setupmatrix[1, j] = (j-1)*gap
			end
			for i in 1:length(s1)+1
				setupmatrix[i,1] = (i-1)*gap
			end
			return setupmatrix
		end

		function nwscorematrix(seq1, seq2; match=1, mismatch=-1, gap=-1)
			scoremat = nwsetupmatrix(seq1, seq2; gap=gap)
			for i in 2:size(scoremat, 1) # iterate through row indices
				for j in 2:size(scoremat, 2) # iterate through column indices
					above = nwscore(seq1[i-1], nothing; match, mismatch, gap) + 							scoremat[i-1,j]
					left = nwscore(nothing, seq2[j-1]; match, mismatch, gap) + 								scoremat[i,j-1]
					diagonal = nwscore(seq1[i-1], seq2[j-1]; match, mismatch, gap) + 						scoremat[i-1,j-1]
					scoremat[i,j] = max(above, left, diagonal)
				end
			end
			return scoremat
		end
end

# ╔═╡ 3bc28dec-531f-4550-8379-2a502d64aed4
function SW(seq1, seq2)
	#use SW alignment function from lab
end

# ╔═╡ 4fe25d29-13a0-458e-98a7-3c7541954b44
function compare(NW,SW)
	#take the Tuple from NW and calculate number of gaps, As, Ts, Cs, and Gs in NW vs SW
	#calculate gc content for each 
	#print comparisons as Tuple 
end

# ╔═╡ 7f495726-3311-4760-8e9a-0d8b618e3927
md"Analysis #2: Convert genome to protein sequences"

# ╔═╡ 27f07737-711a-4ef7-9326-8020781a91e9
function protein_dict()
	#make dictionary to convert DNA --> amino acid sequences 
end

# ╔═╡ 7516a030-b5bb-4b49-ba21-cb93f4f1e24f
function genome_to_protein()
	for seq in 3:genome
		#convert to protein
	end 
end

# ╔═╡ Cell order:
# ╟─9cda8166-e37f-11eb-2e1b-bf74d988c460
# ╟─e7ad56ce-336e-4c1b-9aae-c3a989f7174b
# ╠═a8a8631b-d883-480a-83e6-3015bffb2e3b
# ╠═3bc28dec-531f-4550-8379-2a502d64aed4
# ╠═4fe25d29-13a0-458e-98a7-3c7541954b44
# ╟─7f495726-3311-4760-8e9a-0d8b618e3927
# ╠═27f07737-711a-4ef7-9326-8020781a91e9
# ╠═7516a030-b5bb-4b49-ba21-cb93f4f1e24f
