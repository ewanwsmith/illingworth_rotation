# take SARS-CoV_2 .fasta files, extract all possible open reading frames, and match these to NCBI reference ORFs

#input paths to folders with .fasta files
folder_list = ["data/CAMP000427", "data/CAMP001339", "data/CAMP001490",
                "data/CAMP001523", "data/CAMP002274", "data/CAMP003468",
                "data/CAMP004884", "data/CAMP007136", "data/Kemp"]


#input path to reference ORFs
reference_path = "data/reference/coding_sequences.fasta" 


#load dependencies
using DataFrames
using CSV
using FASTX
using BioSequences


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
        Reference_orf_name = orf_names,
        Sequence = sequences
    )
    df = hcat(df,
    DataFrame(reduce(vcat, permutedims.(split.(df.Reference_orf_name, '|'))),
    [:Gene_accession_and_position, :ORF_name]))
    select!(df, Not(:Reference_orf_name))
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


#find ORFs within consensus sequence
function find_consensus_orfs(consensus_record::FASTX.FASTA.Record, orf_df::DataFrame)
    consensus_sequence = FASTX.sequence(consensus_record)
    
    best_matches = DataFrame(ORF_name=String[], Start_Position=Int[], End_Position=Int[], Matched_Sequence=String[], Match_Length=Int[], Reference_Sequence=String[], Reference_Length=Int[])
    
    for i in 1:nrow(orf_df)
        orf_name = orf_df[i, :ORF_name]
        orf_sequence = orf_df[i, :Sequence]
        
        start_codon = "ATG"
        stop_codons = ["TAA", "TAG", "TGA"]
        
        # Iterate through the consensus sequence to find potential ORFs
        for start in 1:length(consensus_sequence) - length(orf_sequence) + 1
            match_end = start + length(orf_sequence) - 1
            match_sequence = consensus_sequence[start:match_end]
            
            # Check if the potential ORF starts with a start codon
            if startswith(uppercase(match_sequence), start_codon)
                # Check if it ends with a stop codon
                if any(endswith(uppercase(match_sequence), stop_codon) for stop_codon in stop_codons)
                    push!(best_matches, (orf_name, start, match_end, match_sequence, length(match_sequence), orf_sequence, length(orf_sequence)))
                end
            end
        end
    end
    
    return best_matches
end


# match found ORFs to reference ORFs by length and fewest mismatches
function count_mismatches(seq1::AbstractString, seq2::AbstractString) #small function to count mismatches
    return sum(seq1[i] != seq2[i] for i in eachindex(seq1))
end

function match_consensus_orfs(result_df::DataFrame)
    # Calculate n_mismatches for each row
    result_df.n_mismatches = [count_mismatches(row.Matched_Sequence, row.Reference_Sequence) for row in eachrow(result_df)]

    # Group by ORF_name, sort each group by n_mismatches, and take the first row of each group
    best_matches = combine(groupby(result_df, :ORF_name)) do group
        sort!(group, :n_mismatches)
        first(group, 1)
    end

    return best_matches
end

# run above functions
ref_orfs = ref_readin(reference_path)

function find_match_consensus_orfs(folder::String)
    println("Running readin_consensus() function on folder: $folder")
    df = readin_consensus(folder) #readin fasta file

    println("Running find_consensus_orfs() function on folder: $folder")
    orfs_result = find_consensus_orfs(df, ref_orfs) #find possible ORFs in consensus fasta file
    println("Running match_consensus_orfs() function on folder: $folder")
    match_orfs_result = match_consensus_orfs(orfs_result) #find best matches to reference ORFs by length

    csv_path = joinpath(folder, "consensus_ORFs.csv") 
    println("Writing matched_ORFs.csv to folder: $folder")
    CSV.write(csv_path, match_orfs_result) #write results to a CSV in sample folder

    return match_orfs_result
end

#for folder in folder_list
#    find_match_consensus_orfs(folder)
#end

# read variants, find the ORFs they sit in
function locate_variants(folder_path::String)
    # Read Variant_list.csv into variants_df
    variants_path = joinpath(folder_path, "Variant_list.csv")
    variants_df = CSV.read(variants_path, DataFrame)

    # Read Consensus_ORFs.csv into orfs_df
    orfs_path = joinpath(folder_path, "Consensus_ORFs.csv")
    orfs_df = CSV.read(orfs_path, DataFrame)

    # Initialize variant_locations_df
    variant_locations_df = DataFrame(
        ORF_name = String[],
        Start_Position = Int[],
        End_Position = Int[],
        Sequence = String[],
        Variant_Position = Int[],
        Original_Base = String[],
        Variant_Base = String[]
    )

    # Count of no matches
    no_match_count = 0

    # Iterate through each row in variants_df
    for i in 1:size(variants_df, 1)
        position_value = variants_df[i, :Position]

        # Find the row in orfs_df where Position is between Start_Position and End_Position
        matching_row = filter(row -> row.Start_Position <= position_value <= row.End_Position, orfs_df)

        # If a match is found, add a row to variant_locations_df
        if !isempty(matching_row)
            push!(variant_locations_df, (
                matching_row[1, :ORF_name],
                matching_row[1, :Start_Position],
                matching_row[1, :End_Position],
                matching_row[1, :Matched_Sequence],
                variants_df[i, :Position],
                variants_df[i, :Original_Base],
                variants_df[i, :Variant_Base]
            ))
        else
            # If no match is found, increment no_match_count
            no_match_count += 1
        end
    end

    # Print the count of rows with no matches
    println("Number of variants with no matches: $no_match_count")

    return variant_locations_df
end

# pull out codon containing Variant_Base position
function find_codons(df::DataFrame)
    # Create a new column 'Original_Codon'
    df.Original_Codon .= ""

    for i in 1:nrow(df)
        # Check if 'Sequence' column is missing
        if ismissing(df[i, :Sequence])
            continue
        end

        # Extract the sequence, variant position, start position, and original base
        sequence = string(df[i, :Sequence])
        variant_position = df[i, :Variant_Position]
        start_position = df[i, :Start_Position]
        original_base = df[i, :Original_Base]

        # Calculate the adjusted variant position
        adjusted_variant_position = variant_position - start_position

        # Calculate the start position of the selected codon
        codon_start_position = 3 * div(adjusted_variant_position, 3)

        # Extract the codon containing the nth base
        original_codon = sequence[codon_start_position + 1:codon_start_position + 3]

        # Update the 'Original_Codon' column
        df[i, :Original_Codon] = original_codon

        # Check if Original_Base is contained within Original_Codon
        if occursin(original_base, original_codon)
            # Do something if it's contained
        else
            println("Warning: Original_Base is not contained within Original_Codon at row $i.")
        end
    end

    return df
end

kemp_located = locate_variants("data/Kemp")
kemp_codons = find_codons(kemp_located)