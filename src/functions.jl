# functions

#load dependencies
using DataFrames
using CSV
using FASTX
using BioSequences
using DelimitedFiles


# read in reference ORF fasta
function ref_readin(file_path::AbstractString)
    orf_names = String[]
    sequences = String[]
    current_sequence = ""

    open(file_path) do file
        current_orf_name = ""
        for line in eachline(file)
            if occursin(">", line)
                if current_orf_name != ""
                    push!(orf_names, current_orf_name)
                    push!(sequences, current_sequence)
                    current_sequence = ""
                end
                current_orf_name = split(strip(line), '>')[2]
            else
                current_sequence *= strip(line)
            end
        end
        if current_orf_name != ""
            push!(orf_names, current_orf_name)
            push!(sequences, current_sequence)
        end
    end

    df = DataFrame(
        Reference_orf_name=orf_names,
        Sequence=sequences
    )
    df = hcat(df,
        DataFrame(reduce(vcat, permutedims.(split.(df.Reference_orf_name, '|'))),
            [:Gene_accession_and_position, :ORF_name_species]))
    df = hcat(df,
        DataFrame(reduce(vcat, permutedims.(split.(df.ORF_name_species, '['))),
            [:ORF_name, :Species_name]))
    select!(df, Not(:Species_name))
    df = hcat(df,
        DataFrame(reduce(vcat, permutedims.(split.(df.Gene_accession_and_position, ':'))),
            [:Accession, :Position]))
    select!(df, Not(:Gene_accession_and_position))
    df = hcat(df,
        DataFrame(reduce(vcat, permutedims.(split.(df.Position, ".."))),
            [:Start_Position, :End_Position]))
    select!(df, Not(:Position))
    select!(df, [:ORF_name, :Start_Position, :End_Position, :Sequence])
    df.Start_Position = parse.(Int64, df.Start_Position)
    df.End_Position = parse.(Int64, df.End_Position)
    return df
end


# read in consensus sequences
function readin_consensus(folder::AbstractString)
    # Construct the full path to the Consensus.fa file
    file_path = joinpath(folder, "Consensus.fa")

    # Check if the file exists
    if isfile(file_path)
        # Read the FASTA file using FASTX.FASTA.Reader
        consensus_records = FASTX.FASTA.Reader(open(file_path))

        # Check if there's at least one record
        if !eof(consensus_records)
            # Extract the first record
            consensus_sequence = read(consensus_records)

            close(consensus_records)  # Close the file handle

            return consensus_sequence
        else
            error("Consensus.fa file is empty.")
        end
    else
        error("Consensus.fa file not found in the specified folder.")
    end
end

# pull protein sequences from positions
function pull_proteins(record::FASTX.FASTA.Record, positions_df::DataFrame)
    # Get sequence from the FASTA record
    sequence = FASTX.sequence(record)
    
    # Initialize an empty DataFrame to store the results
    result_df = DataFrame(
        Protein = String[],
        Start = Int[],
        End = Int[],
        Sequence = String[]
    )
    
    # Iterate through each row in the positions DataFrame
    for row in eachrow(positions_df)
        protein_name = row.Protein
        position_str = row.Position
        
        # Extract start and end positions from the position string
        start_end_pairs = split(position_str, ',')
        
        for pair in start_end_pairs
            # Convert the substring to a string explicitly
            pair_string = String(pair)

            match_positions = match(r"(\d+)..(\d+)", pair_string)
            start_position = parse(Int, match_positions.captures[1])
            end_position = parse(Int, match_positions.captures[2])

            # Extract the subsequence from start_position to end_position
            subsequence = sequence[start_position:end_position]

            # Append the result to the DataFrame
            push!(result_df, (protein_name, start_position, end_position, subsequence))
        end
    end
    
    return result_df
end

# join together frames of genes with frameshifts
function join_frames(df)
    # Sort the DataFrame by Protein, Start, and End
    sort!(df, [:Protein, :Start, :End])

    # Initialize variables to store the joined rows
    joined_proteins = String[]
    joined_start = nothing
    joined_end = nothing
    joined_sequence = ""

    # Iterate through the DataFrame to join rows with the same Protein
    for row in eachrow(df)
        if !isempty(joined_proteins) && row.Protein == joined_proteins[end]
            # If the current row has the same Protein as the previous one, update End and concatenate Sequence
            joined_end = max(joined_end, row.End)
            joined_sequence *= row.Sequence
        else
            # If the current row has a different Protein, store the joined row and reset variables
            if !isempty(joined_proteins)
                df[df.Protein .== joined_proteins[end], :Start] .= joined_start
                df[df.Protein .== joined_proteins[end], :End] .= joined_end
                df[df.Protein .== joined_proteins[end], :Sequence] .= joined_sequence
            end

            push!(joined_proteins, row.Protein)
            joined_start = row.Start
            joined_end = row.End
            joined_sequence = row.Sequence
        end
    end

    # Update the last joined row after the loop
    if !isempty(joined_proteins)
        df[df.Protein .== joined_proteins[end], :Start] .= joined_start
        df[df.Protein .== joined_proteins[end], :End] .= joined_end
        df[df.Protein .== joined_proteins[end], :Sequence] .= joined_sequence
    end

    # Drop duplicate rows (keeping the first occurrence)
    df = unique(df, [:Protein])
end

# run above functions
function find_genes(folder::String)
    println("Running readin_consensus() function on folder: $folder")
    df = readin_consensus(folder) #readin fasta file

    println("Running pull_proteins() function on folder: $folder")
    df = pull_proteins(df, positions_df) #pull genes out from consensus sequences
    println("Running join_frames() function on folder: $folder")
    df = join_frames(df) #join frames across frameshifts

    csv_path = joinpath(folder, "pulled_genes.csv") 
    println("Writing pulled_genes.csv to folder: $folder")
    CSV.write(csv_path, df) #write results to a CSV in sample folder

    return df
end

function readin_variants(folder_path::AbstractString)
    # Construct the full path to the Variant_list.csv file
    file_path = joinpath(folder_path, "Variant_list.csv")
    
    try
        # Try to read the CSV file into a DataFrame
        variants_df = CSV.File(file_path) |> DataFrame
        
        # Convert Variant_Base and Original_Base to String1
        variants_df[!, :Variant_Base] .= string.(variants_df[!, :Variant_Base])
        variants_df[!, :Original_Base] .= string.(variants_df[!, :Original_Base])
        
        # Replace "true" with "T" in Original_Base and Variant_Base columns
        variants_df[!, :Original_Base] .= replace.(variants_df[!, :Original_Base], "true" => "T")
        variants_df[!, :Variant_Base] .= replace.(variants_df[!, :Variant_Base], "true" => "T")
        
        return variants_df
    catch e
        # Handle the case when the file is not found or there is an error in reading
        println("Error: ", e)
        println("Could not read Variant_list.csv. Please check the file path.")
        return DataFrame()  # Return an empty DataFrame in case of an error
    end
end


# read variants, find the ORFs they sit in
function locate_variants(variants_df::DataFrame, frames_df::DataFrame)
    # Initialize variant_locations_df
    variant_locations_df = DataFrame(
        Protein = String[],
        Start_Position = Int[],
        End_Position = Int[],
        Sequence = String[],
        Variant_Position = Int[],
        Original_Base = String[],
        Variant_Base = String[]
    )

    # Count of no matches
    no_match_count = 0

    # Sort frames_df by Start column
    sorted_frames_df = sort(frames_df, [:Protein, :Start])

    # Iterate through each row in variants_df
    for i in 1:size(variants_df, 1)
        position_value = variants_df[i, :Position]

        # Find the rows in sorted_frames_df where Position is between Start and End
        matching_rows = filter(row -> row.Start <= position_value <= row.End, sorted_frames_df)

        # If matches are found, join them
        if !isempty(matching_rows)
            protein_value = matching_rows[1, :Protein]
            start_value = matching_rows[1, :Start]
            end_value = matching_rows[end, :End]  # Take the last numerical End value
            sequence_value = join(matching_rows[!, :Sequence], "")  # Join Sequence values

            variant_base = string(variants_df[i, :Variant_Base])

            push!(variant_locations_df, (
                protein_value,
                start_value,
                end_value,
                sequence_value,
                (variants_df[i, :Position] + 2),
                variants_df[i, :Original_Base],
                variant_base
            ))
        else
            no_match_count += 1
        end
    end

    println("Number of variants with no matches: $no_match_count")

    return variant_locations_df
end

# create variant_sequence
function substitute_variants(dataframe::DataFrame)
    # Create a new column for Base_Position
    dataframe.Adj_Variant_Position = dataframe.Variant_Position .- dataframe.Start_Position

    # Create a new column for Variant_Sequence
    dataframe.Variant_Sequence = Vector{String}(undef, nrow(dataframe))

    # Iterate over each row in the DataFrame
    for i in 1:nrow(dataframe)
        # Extract the original sequence
        sequence = string(dataframe[i, :Sequence])

        # Extract the variant position
        adjusted_variant_position = dataframe[i, :Adj_Variant_Position]

        # Check if the variant position is within the sequence length
        if 1 <= adjusted_variant_position <= length(sequence)
            # Extract the variant base
            variant_base = string(dataframe[i, :Variant_Base])

            # Check if the character at Base_Position is equal to Original_Base
            original_base = string(sequence[adjusted_variant_position])
            if original_base != dataframe[i, :Original_Base]
                println("Warning: Original_Base in row $i does not match base at position $adjusted_variant_position in the Sequence.")
            end

            # Create the new sequence with the substitution
            new_sequence = string(sequence[1:adjusted_variant_position-1], variant_base, sequence[adjusted_variant_position+1:end])

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
        sequence = string(df[row, :Sequence])
        base_position = df[row, :Adj_Variant_Position]
        
        # Check if base_position is valid
        if base_position < 1 || base_position > length(sequence)
            throw(ArgumentError("Invalid Base_Position for row $row"))
        end
        
        # Extract the codon based on the Adj_Variant_Position
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
        base_position = df[row, :Adj_Variant_Position]
        
        # Check if base_position is valid
        if base_position < 1 || base_position > length(sequence)
            throw(ArgumentError("Invalid Base_Position for row $row"))
        end
        
        # Extract the codon based on the Adj_Variant_Position
        codon_index = (base_position - 1) ÷ 3 + 1
        push!(codons, split_into_codons(sequence)[codon_index])
    end
    
    # Add the codons as a new column in the DataFrame
    df[!, :Variant_Codon] = codons
    
    return df
end

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

#translate original & variant codons & determine synonymity
function translate_codons(df)
    # Translate Original_Codon and Variant_Codon to amino acids
    df.Original_AA = map(codon_to_aa, df.Original_Codon)
    df.Variant_AA = map(codon_to_aa, df.Variant_Codon)

    # Determine if the amino acids are synonymous
    df.Is_Synonymous = ifelse.(df.Original_AA .== df.Variant_AA, "Yes", "No")

    return df
end

#T bases keep being mistaken for "True", so explicitly correct
function true_to_T(df::DataFrame)
    # Replace "true" with "T" in Original_Base column
    df[!, :Original_Base] .= replace.(df[!, :Original_Base], "true" => "T")
    
    # Replace "true" with "T" in Variant_Base column
    df[!, :Variant_Base] .= replace.(df[!, :Variant_Base], "true" => "T")
    
    return df
end

# run above functions
function pull_translate_codons(folder::String)
    println("Running readin_consensus() function on folder: $folder")
    df = readin_consensus(folder) #readin fasta file
    println("Running readin_variants() function on folder: $folder")
    var = readin_variants(folder) #readin variants file
    println("Running pull_proteins() function on folder: $folder")
    df = pull_proteins(df, positions_df) #pull genes out from consensus sequences
    println("Running locate_variants() function on folder: $folder")
    df = locate_variants(var, df) #find variants in ORFs

    println("running true_to_T() function on folder: $folder")
    df = true_to_T(df) # fix true / T issue 

    println("Running substitute_variants() function on folder: $folder")
    df = substitute_variants(df) #create variant_sequence
    println("Running find_original_codons() function on folder: $folder")
    df = find_original_codons(df) #pull original codon
    println("Running find_variant_codons() function on folder: $folder")
    df = find_variant_codons(df) #pull variant codon
    println("Running translate_codons() function on folder $folder")
    df = translate_codons(df)

    csv_path = joinpath(folder, "codons.csv") 
    println("writing codons.csv to folder: $folder")
    CSV.write(csv_path, df) #write results to a CSV in sample folder

    return df
end

function find_rates(folder_path::AbstractString)
    # Construct the full path to the Mean_rates.txt file
    rates_path = joinpath(folder_path, "Mean_rates.txt")
    probs_path = joinpath(folder_path, "Fixation_probabilities.txt")

    try
        # Try reading the files with CSV.File and tab delimiter
        rates_df = CSV.File(rates_path, delim='\t', types=[Int, String, String, Float64], header=["Position", "Original_Base", "Variant_Base", "Evo_rate"]) |> DataFrame
        probs_df = CSV.File(probs_path, delim='\t', types=[Int, String, String, Float64], header=["Position", "Original_Base", "Variant_Base", "Pr_fixation"]) |> DataFrame

        return rates_df, probs_df
    catch e
        # If there is an error, print an error message
        println("Error reading the files: $e")
    end

    # If no successful read, return empty DataFrames
    return DataFrame(), DataFrame()
end


function join_rates(rates_df::DataFrame, probs_df::DataFrame)
    try
        # Join DataFrames based on the "Position" column using inner join
        merged_df = innerjoin(rates_df, probs_df, on=:Position, makeunique=true)

        # Drop the unwanted columns
        select!(merged_df, Not(:Original_Base_1, :Variant_Base_1))
        
        return merged_df
    catch e
        println("Error joining DataFrames: $e")
        return DataFrame()  # Return an empty DataFrame in case of an error
    end
end

# add evolution rates and Pr(fixation) data from model
function add_new_data(folder_path::AbstractString)
    # Call find_rates to get rates_df and probs_df
    rates_df, probs_df = find_rates(folder_path)

    # Check if find_rates encountered an error
    if isempty(rates_df) || isempty(probs_df)
        return DataFrame()  # Return an empty DataFrame in case of an error
    end

    try
        # Join DataFrames based on the "Position" column using inner join
        merged_df = innerjoin(rates_df, probs_df, on=:Position, makeunique=true)

        # Drop the unwanted columns
        select!(merged_df, Not(:Original_Base_1, :Variant_Base_1))

        # Construct the full path to the codons.csv file
        codons_path = joinpath(folder_path, "codons.csv")

        # Read the codons.csv file into a DataFrame
        codons_df = CSV.File(codons_path) |> DataFrame

        # Adjust the Position values in the merged_df
        merged_df.Position .= merged_df.Position .+ 2

        # Join DataFrames based on the adjusted "Position" and "Variant_Position" columns
        final_merged_df = innerjoin(merged_df, codons_df, on=:Position => :Variant_Position, makeunique=true)

        # Drop the additional columns
        select!(final_merged_df, Not(:Original_Base_1, :Variant_Base_1))

        # Display the modified DataFrame
        display(final_merged_df)

        # Save the final merged DataFrame as codons.csv in the specified folder
        CSV.write(codons_path, final_merged_df)

        return final_merged_df
    catch e
        println("Error reading or joining DataFrames: $e")
        return DataFrame()  # Return an empty DataFrame in case of an error
    end
end
