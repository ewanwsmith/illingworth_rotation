# take SARS-CoV_2 .fasta files, extract all possible open reading frames, and match these to NCBI reference ORFs

#input paths to folders with .fasta files
folder_list = ["data/CAMP000427", "data/CAMP001339", "data/CAMP001490",
                "data/CAMP001523", "data/CAMP002274", "data/CAMP003468",
                "data/CAMP004884", "data/CAMP007136", "data/Kemp"]


#input path to reference ORFs
reference_path = "data/reference/coding_sequences.fasta" 

#load functions
include("/Users/ewansmith/Documents/PhD/Rotation 1 - Illingworth/illingworth_rotation/src/functions.jl")

# readin reference ORFs
ref_orfs = ref_readin(reference_path)

# find ORFs within consensus sequences
for folder in folder_path
    find_match_consensus_orfs(folder)
end

# locate variants within ORFs, pull out affected codons & deterine synonymity
for folder in folder_path
    pull_translate_codons(folder)
end