# take SARS-CoV_2 .fasta files, extract all possible open reading frames, and match these to NCBI reference ORFs

# input paths to folders with .fasta files
folder_list = ["data/CAMP000427", "data/CAMP001339", "data/CAMP001490",
                "data/CAMP001523", "data/CAMP002274", "data/CAMP003468",
                "data/CAMP004884", "data/CAMP007136", "data/Kemp"]


# input path to reference files
reference_orf_path = "data/reference/coding_sequences.fasta" 
reference_protein_path = "data/reference/protein_sequences.fasta"
positions_path = "data/reference/protein_positions.csv"

# load functions
include("/Users/ewansmith/Documents/PhD/Rotation 1 - Illingworth/illingworth_rotation/src/functions.jl")

# readin reference ORFs
ref_orfs = ref_readin(reference_orf_path)

# readin protein positions
positions_df = CSV.File(positions_path) |> DataFrame

# find ORFs within consensus sequences
for folder in folder_list
    find_genes(folder)
end

# locate variants within ORFs, pull out affected codons & deterine synonymity
for folder in folder_list
    pull_translate_codons(folder)
end

# add model outputs to codons.csv, save as new codons.csv
#for folder in folder_list
#    add_new_data(folder)
#end