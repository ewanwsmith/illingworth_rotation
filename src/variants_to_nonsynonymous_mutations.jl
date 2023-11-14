# input paths to folders with Variant_list.out files
folder_list = ["data/CAMP000427", "data/CAMP001339", "data/CAMP001490",
"data/CAMP001523", "data/CAMP002274", "data/CAMP003468",
"data/CAMP004884", "data/CAMP007136", "data/Kemp"] 

# input path to reference ORFs
reference_path = "data/reference/coding_sequences.csv"

# load dependencies
using DataFrames
using CSV

# find which ORF has each variant
function locate_variants(folder::AbstractString, reference_orf_path::AbstractString)
    # Read Variant_list.csv into a DataFrame
    variant_df = CSV.File(joinpath(folder, "Variant_list.csv"), types=[Int64, Char, Char], header=true) |> DataFrame

    # Read reference ORFs CSV into a DataFrame
    reference_orf_df = CSV.File(reference_orf_path) |> DataFrame

    # Initialize an empty DataFrame to store results
    result_df = DataFrame(Variant_Position = Int[], Original_Base = Char[], Variant_Base = Char[],
                          Start_Position = Int[], End_Position = Int[], ORF_Name = String[], Sequence = String[])

    # Iterate over each variant in Variant_list.csv
    for row in eachrow(variant_df)
        variant_position = row.Position
        original_base = row.Original_Base
        variant_base = row.Variant_Base

        # Iterate over each ORF in reference ORFs DataFrame
        for orf_row in eachrow(reference_orf_df)
            start_position = orf_row.Start_Position
            end_position = orf_row.End_Position
            orf_name = orf_row.ORF_name
            sequence = orf_row.Sequence

            # Check if the variant_position is within the ORF range
            if start_position <= variant_position <= end_position
                push!(result_df, (variant_position, original_base, variant_base, start_position, end_position, orf_name, sequence))
            end
        end
    end

    return result_df
end

using DataFrames

function substitute_variants(dataframe::DataFrame)
    # Create a new column for Variant_Sequence
    dataframe.Variant_Sequence = Vector{String}(undef, nrow(dataframe))

    # Iterate over each row in the DataFrame
    for i in 1:nrow(dataframe)
        # Calculate the modified variant position
        variant_position = dataframe[i, :Variant_Position] - dataframe[i, :Start_Position] + 1
        variant_base = string(dataframe[i, :Variant_Base])

        # Extract the original sequence
        sequence = string(dataframe[i, :Sequence])

        # Check if the variant position is within the sequence length
        if 1 <= variant_position <= length(sequence)
            # Create the new sequence with the substitution
            new_sequence = string(sequence[1:variant_position-1], variant_base, sequence[variant_position+1:end])

            # Update the Variant_Sequence column
            dataframe[i, :Variant_Sequence] = new_sequence
        end
    end

    return dataframe
end

# find variant codons
using DataFrames

function find_codons(df::DataFrame)
    # Calculate the base position based on Variant_Position and Start_Position
    df.Base_Position .= df.Variant_Position .- df.Start_Position
    
    # Function to extract the three-base codon based on counting in threes
    function extract_codon(seq, position)
        codon_start = (position - 1) ÷ 3 * 3 + 1
        codon_end = min(length(seq), codon_start + 2)
        return seq[codon_start:codon_end]
    end
    
    # Extract Original_Base and Variant_Base using the extract_codon function
    df[!, :Original_Codon] .= extract_codon.(df[:, :Sequence], df[:, :Base_Position])
    df[!, :Variant_Codon] .= extract_codon.(df[:, :Variant_Sequence], df[:, :Base_Position])
    
    return df
end

function translate_codons(df)
    # Function to map codons to amino acids
    function codon_to_aa(codon)
        codon_dict = Dict("TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
                          "CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
                          "ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
                          "GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
                          "TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
                          "CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
                          "ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
                          "GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
                          "TAT" => "Y", "TAC" => "Y", "TAA" => "*", "TAG" => "*",
                          "CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
                          "AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
                          "GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
                          "TGT" => "C", "TGC" => "C", "TGA" => "*", "TGG" => "W",
                          "CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
                          "AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
                          "GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G")
        
        return get(codon_dict, codon, "Unknown")
    end

    # Translate Original_Codon and Variant_Codon to amino acids
    df.Original_AA = map(codon_to_aa, df.Original_Codon)
    df.Variant_AA = map(codon_to_aa, df.Variant_Codon)

    # Determine if the amino acids are synonymous
    df.Is_Synonymous = ifelse.(df.Original_AA .== df.Variant_AA, "Yes", "No")

    return df
end

# run above functions
for folder in folder_list
    println("running locate_variants() function on folder: $folder")
    variant_locations = locate_variants(folder, reference_path)
    println("running substitute_variants() function on folder: $folder")
    variant_sequences = substitute_variants(variant_locations)
    println("running find_codons() function on folder: $folder")
    variant_codons = find_codons(variant_sequences)
    println("running translate_codons() function on folder: $folder")
    variant_translated = translate_codons(variant_codons)

    println("saving dataframe as variant_sequences.csv in folder: $folder")
    CSV.write(joinpath(folder, "variant_sequences.csv"), variant_translated)
end
