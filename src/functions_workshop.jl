# somewhere to put functions as they're worked on

#load functions
include("/Users/ewansmith/Documents/PhD/Rotation 1 - Illingworth/illingworth_rotation/src/functions.jl")

# read in protein positions
positions_path = ("data/reference/protein_positions.csv")

positions_df = CSV.File(positions_path) |> DataFrame


# get Kemp consensus
kemp = readin_consensus("data/Kemp")

# pull protein sequences from positions
using DataFrames
using FASTX

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

kemp = pull_proteins(kemp, positions_df)

function join_frames(df)
    # Assuming your data is stored in a DataFrame called df

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

kemp = join_frames(kemp)