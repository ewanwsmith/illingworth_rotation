# somewhere to put functions as they're worked on

using DataFrames
using Dates
using DelimitedFiles
using CSV

# translate dictionary ORFs
function translate_orfs(fasta_dict::Dict{String, Vector{String}})
    translated_orfs = Dict{String, Vector{String}}()

    codon_table = Dict(
        "TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
        "TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
        "TAT" => "Y", "TAC" => "Y", "TAA" => "*", "TAG" => "*",
        "TGT" => "C", "TGC" => "C", "TGA" => "*", "TGG" => "W",
        "CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
        "CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
        "CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
        "CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
        "ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
        "ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
        "AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
        "AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
        "GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
        "GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
        "GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
        "GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G"
    )

    for (header, sequences) in fasta_dict
        translated_sequences = Vector{String}()
        for seq in sequences
            codons = [seq[i:i+2] for i in 1:3:length(seq)-2]
            translated_sequence = join([get(codon_table, codon, "X") for codon in codons], "")
            push!(translated_sequences, translated_sequence)
        end
        translated_orfs[header] = translated_sequences
    end

    return translated_orfs
end

matched427 = CSV.read("data/CAMP000427/matched_orfs.csv", DataFrame)
matchedkemp = CSV.read("data/kemp/matched_orfs.csv", DataFrame)


# pullout mismatches between first sequence in time / Px and reference sequence

function find_first_variants(df::DataFrame)
    result_df = DataFrame(
        Sequence_Name = String[],
        ORF_name = String[],
        Original_Base = String[],
        Variant_Base = String[]
    )

    # Group the DataFrame by ORF_name
    grouped_df = groupby(df, :ORF_name)

    for sub_df in grouped_df
        # Sort the sub-dataframe by Date
        sub_df = sort(sub_df, order(:Date))

        # Extract the first Matched Sequence
        first_matched_sequence = sub_df.Matched_Sequence[1]

        reference_sequence = sub_df.Reference_Sequence[1]
        start_position = sub_df.Start_Position[1]

        for (i, (ref_base, match_base)) in enumerate(zip(reference_sequence, first_matched_sequence))
            if ref_base != match_base
                original_base_position = "$(start_position + i):$ref_base"
                variant_base_position = "$(start_position + i):$match_base"
                push!(result_df, (sub_df.Sequence_Name[1], sub_df.ORF_name[1], original_base_position, variant_base_position))
            end
        end
    end

    return result_df
end

firstvariants427 = find_first_variants(matched427)
firstvariantskemp = find_first_variants(matchedkemp)

#pullout variants in 2:end samples 
function find_subs_variants(df::DataFrame)
    result_df = DataFrame(
        Sequence_Name = String[],
        ORF_name = String[],
        Population = String[],
        Original_Base = String[],
        Variant_Base = String[]
    )

    # Group the DataFrame by ORF_name and Population
    grouped_df = groupby(df, [:ORF_name, :Population])

    for sub_df in grouped_df
        # Sort the sub-dataframe by Date
        sub_df = sort(sub_df, order(:Date))

        for i in 2:size(sub_df, 1)
            previous_sequence = sub_df.Matched_Sequence[i - 1]
            current_sequence = sub_df.Matched_Sequence[i]
            start_position = sub_df.Start_Position[i]

            for (j, (prev_base, curr_base)) in enumerate(zip(previous_sequence, current_sequence))
                if prev_base != curr_base
                    original_base_position = "$(start_position + j):$prev_base"
                    variant_base_position = "$(start_position + j):$curr_base"
                    push!(result_df, (sub_df.Sequence_Name[i], sub_df.ORF_name[i], sub_df.Population[i], original_base_position, variant_base_position))
                end
            end
        end
    end

    return result_df
end

subsvariants427 = find_subs_variants(matched427)
subsvariantskemp = find_subs_variants(matchedkemp)
