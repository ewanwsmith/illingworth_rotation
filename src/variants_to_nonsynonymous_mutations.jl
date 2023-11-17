#load dependencies
include("/Users/ewansmith/Documents/PhD/Rotation 1 - Illingworth/illingworth_rotation/src/fasta_to_matched_ORFs.jl")

# build consensus sequences for each ORF, save as .csv in folder
function build_consensus(folder_path::AbstractString)
    # Load the CSV file
    csv_path = joinpath(folder_path, "matched_orfs.csv")
    df = CSV.File(csv_path) |> DataFrame
    
    # Group by ORF_name
    grouped_df = groupby(df, :ORF_name)
    
    # Initialize an empty DataFrame to store consensus sequences
    consensus_df = DataFrame(ORF_name = String[],
                             Start_Position = Int[], 
                             End_Position = Int[], 
                             Reference_Sequence = String[],
                             Consensus_Sequence = String[])
    
    # Iterate over each group
    for (group_idx, group) in enumerate(grouped_df)
        orf_name = first(group[!, :ORF_name])
        
        # Sort the group by Date
        sorted_group = sort(group, order(:Date))
        
        # Initialize consensus sequence with the first sequence
        consensus_seq = collect(first(sorted_group.Matched_Sequence))

        # Iterate over each position in the sequences
        for i in 2:length(sorted_group.Matched_Sequence)
            # Check if the current position has an 'N'
            if consensus_seq[i] == 'N'
                # Iterate over the sequences to find a non-'N' base
                for row in eachrow(sorted_group)
                    if row.Matched_Sequence[i] != 'N'
                        consensus_seq[i] = row.Matched_Sequence[i]
                        break
                    end
                end
            end
        end

        # Convert the consensus sequence back to a string
        consensus_seq_str = join(consensus_seq)
        
        # Retrieve other relevant information
        start_pos = first(sorted_group.Start_Position)
        end_pos = first(sorted_group.End_Position)
        ref_seq = first(sorted_group.Reference_Sequence)
        
        # Add the consensus sequence to the DataFrame
        push!(consensus_df, (ORF_name = orf_name, 
                             Consensus_Sequence = consensus_seq_str,
                             Start_Position = start_pos, 
                             End_Position = end_pos, 
                             Reference_Sequence = ref_seq))
        
    end
    
    # Save the consensus DataFrame to a CSV file
    consensus_path = joinpath(folder_path, "Consensus_ORFs.csv")
    CSV.write(consensus_path, consensus_df)
    
    return consensus_df
end


# find which ORF has each variant
function locate_variants(folder_path::AbstractString, consensus_df::DataFrame)
    # Read Variant_list.csv into a DataFrame
    variant_list_path = joinpath(folder_path, "Variant_list.csv")
    variant_df = CSV.File(variant_list_path, types=Dict(:Position => Int, :Original_Base => Char, :Variant_Base => Char)) |> DataFrame

    # Initialize an empty DataFrame to store the results
    result_df = DataFrame(
        Variant_Position = Int[],
        Original_Base = Char[],
        Variant_Base = Char[],
        Start_Position = Int[],
        End_Position = Int[],
        ORF_Name = String[],
        Consensus_Sequence = String[]
    )

    # Iterate over each row in the Variant_list DataFrame
    for row_v in eachrow(variant_df)
        variant_position = row_v.Position

        # Iterate over each row in the consensus_df DataFrame
        for row_o in eachrow(consensus_df)
            # Check if the variant_position is within the range of the open reading frame
            if variant_position >= row_o.Start_Position && variant_position <= row_o.End_Position
                # Extract relevant information and append to result_df
                push!(result_df, (
                    variant_position,
                    row_v.Original_Base[1],
                    row_v.Variant_Base[1],
                    row_o.Start_Position,
                    row_o.End_Position,
                    row_o.ORF_name,
                    row_o.Consensus_Sequence
                ))
                break  # Break the inner loop since we found a match
            end
        end
    end

    return result_df
end


using DataFrames

# create full sequence with variant edited in
function substitute_variants(dataframe::DataFrame)
    # Create a new column for Base_Position
    dataframe.Base_Position = dataframe.Variant_Position .- dataframe.Start_Position

    # Create a new column for Variant_Sequence
    dataframe.Variant_Sequence = Vector{String}(undef, nrow(dataframe))

    # Iterate over each row in the DataFrame
    for i in 1:nrow(dataframe)
        # Extract the original sequence
        sequence = string(dataframe[i, :Consensus_Sequence])

        # Extract the variant position
        variant_position = dataframe[i, :Base_Position]

        # Check if the variant position is within the sequence length
        if 1 <= variant_position <= length(sequence)
            # Extract the variant base
            variant_base = string(dataframe[i, :Variant_Base])

            # Create the new sequence with the substitution
            new_sequence = string(sequence[1:variant_position-1], variant_base, sequence[variant_position+1:end])

            # Update the Variant_Sequence column
            dataframe[i, :Variant_Sequence] = new_sequence
        end
    end

    return dataframe
end

# find original codons
function find_original_codons(df::DataFrame)
    # Function to split a DNA sequence into codons
    function split_into_codons(sequence)
        return [sequence[i:i+2] for i in 1:3:length(sequence)-2]
    end
    
    # Process each row in the DataFrame
    codons = []
    for row in 1:size(df, 1)
        sequence = string(df[row, :Consensus_Sequence])
        base_position = df[row, :Base_Position]
        
        # Check if base_position is valid
        if base_position < 1 || base_position > length(sequence)
            throw(ArgumentError("Invalid Base_Position for row $row"))
        end
        
        # Extract the codon based on the Base_Position
        codon_index = (base_position - 1) ÷ 3 + 1
        push!(codons, split_into_codons(sequence)[codon_index])
    end
    
    # Add the codons as a new column in the DataFrame
    df[!, :Original_Codon] = codons
    
    return df
end

# find variant codons
function find_variant_codons(df::DataFrame)
    # Function to split a DNA sequence into codons
    function split_into_codons(sequence)
        return [sequence[i:i+2] for i in 1:3:length(sequence)-2]
    end
    
    # Process each row in the DataFrame
    codons = []
    for row in 1:size(df, 1)
        sequence = string(df[row, :Variant_Sequence])
        base_position = df[row, :Base_Position]
        
        # Check if base_position is valid
        if base_position < 1 || base_position > length(sequence)
            throw(ArgumentError("Invalid Base_Position for row $row"))
        end
        
        # Extract the codon based on the Base_Position
        codon_index = (base_position - 1) ÷ 3 + 1
        push!(codons, split_into_codons(sequence)[codon_index])
    end
    
    # Add the codons as a new column in the DataFrame
    df[!, :Variant_Codon] = codons
    
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
    println("running build_consensus() function on folder: $folder")
    consensus_sequences = build_consensus(folder)
    println("running locate_variants() function on folder: $folder")
    variant_locations = locate_variants(folder, consensus_sequences)
    println("running substitute_variants() function on folder: $folder")
    variant_sequences = substitute_variants(variant_locations)
    println("running find_original_codons() function on folder: $folder")
    original_codons = find_original_codons(variant_sequences)
    println("running find_variant_codons() function on folder: $folder")
    variant_codons = find_variant_codons(original_codons)
    println("running translate_codons() function on folder: $folder")
    variant_translated = translate_codons(variant_codons)

    println("saving dataframe as variant_sequences.csv in folder: $folder")
    CSV.write(joinpath(folder, "variant_sequences.csv"), variant_translated)
end

kemp_consensus = build_consensus("data/Kemp")
kemp_located = locate_variants("data/Kemp", kemp_consensus)
kemp_substituted = substitute_variants(kemp_located)
kemp_original_codons = find_original_codons(kemp_substituted)
kemp_variant_codons = find_variant_codons(kemp_original_codons)
kemp_translated = translate_codons(kemp_variant_codons)