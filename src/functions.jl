# somewhere to put functions as they're worked on

using CSV
using DataFrames

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

# find codon containing variant
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
