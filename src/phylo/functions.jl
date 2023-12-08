# functions for trying phylogenetics

#input paths to folders with .fasta files
folder_list = ["data/CAMP000427", "data/CAMP001339", "data/CAMP001490",
                "data/CAMP001523", "data/CAMP002274", "data/CAMP003468",
                "data/CAMP004884", "data/CAMP007136", "data/Kemp"]


# load functions
include("/Users/ewansmith/Documents/PhD/Rotation 1 - Illingworth/illingworth_rotation/src/functions.jl")

# pull out consensus fasta files to a single fasta 
using BioSequences
using FASTX

function pull_consensus(folders::Vector{String}, output_file::String)
    consensus_records = BioSequences.FASTX.Record[]

    for folder in folders
        consensus_records_in_folder = readin_consensus(folder)

        # Merge the FASTX records
        append!(consensus_records, consensus_records_in_folder)
    end

    # Write the merged records to a .fa file
    BioSequences.FASTX.write(output_file, consensus_records)
end




pull_consensus(folder_list, "data/phylo/pulled_consensus.fa")