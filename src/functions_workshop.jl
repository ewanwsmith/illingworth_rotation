# somewhere to put functions as they're worked on

#load functions
include("/Users/ewansmith/Documents/PhD/Rotation 1 - Illingworth/illingworth_rotation/src/functions.jl")


using DataFrames

function find_protein(df::DataFrame, fasta_path::AbstractString)
    # Read reference sequences from the fasta file
    reference_sequences = Dict{String, String}()

    # Parse the fasta file
    current_id = ""
    current_sequence = ""
    for line in eachline(fasta_path)
        if startswith(line, '>')
            if !isempty(current_id)
                reference_sequences[current_id] = current_sequence
            end
            current_id = split(strip(line[2:end]), ' ')[1]
            current_sequence = ""
        else
            current_sequence *= strip(line)
        end
    end
    reference_sequences[current_id] = current_sequence

    # Function to find the protein theme
    function get_protein_theme(aa_sequence, variant_aa_sequence)
        variant_position = findfirst(x -> x[1] != x[2], zip(collect(aa_sequence), collect(variant_aa_sequence)))

        if variant_position === nothing
            return "Unknown"
        end

        variant_position = first(variant_position)

        # Find the protein theme
        for (id, ref_seq) in pairs(reference_sequences)
            ref_chars = collect(ref_seq)
            aa_chars = collect(aa_sequence)
            variant_chars = collect(variant_aa_sequence)

            if variant_position <= length(ref_chars) && ref_chars[variant_position] == variant_chars[variant_position]
                return id
            end
        end

        return "Unknown"
    end

    # Apply the function to the DataFrame and add the "Protein" column
    df.Protein = get_protein_theme.(df.AA_Sequence, df.Variant_AA_Sequence)

    # Print the DataFrame with the new "Protein" column
    println(df)

    return df
end


