# somewhere to put functions as they're worked on

# convert Times.in to Dates.in 
using DataFrames
using Dates
using DelimitedFiles
using CSV

function convert_to_date(folder::String)
    times_file = joinpath(folder, "Times.in")
    dates_file = joinpath(folder, "Dates.in")
    if !isfile(times_file)
        error("File 'Times.in' not found in the folder.")
    end

    times = readdlm(times_file, '\t', Int)
    dates = Date[]
    for time in times
        push!(dates, Date(1900, 1, 1) + Day(time))
    end

    df = DataFrame(Dates = dates)
    CSV.write(dates_file, df)

    return df
end

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


using CSV 
using DataFrames

matched427 = CSV.read("data/CAMP000427/matched_orfs.csv", DataFrame)
matchedkemp = CSV.read("data/kemp/matched_orfs.csv", DataFrame)

# try for a variant finding function
function find_variants(df::DataFrame)
    result_df = DataFrame(
        Sequence_Name = String[],
        ORF_name = String[],
        Original_Base_Position = String[],
        Variant_Base_Position = String[]
    )

    for row in eachrow(df)
        sequence_name = row.Sequence_Name
        orf_name = row.ORF_name
        reference_sequence = row.Reference_Sequence
        matched_sequence = row.Matched_Sequence
        start_position = row.Start_Position

        for (i, (ref_base, match_base)) in enumerate(zip(reference_sequence, matched_sequence))
            if ref_base != match_base
                original_base_position = "$(start_position + i):$ref_base"
                variant_base_position = "$(start_position + i):$match_base"
                push!(result_df, (sequence_name, orf_name, original_base_position, variant_base_position))
            end
        end
    end

    return result_df
end

variants427 = find_variants(matched427)
variantskemp = find_variants(matchedkemp)

using DataFrames

function count_variants(df::DataFrame)
    # Count rows per value of ORF_name
    variants_by_ORF = combine(groupby(df, :ORF_name), nrow)

    # Count rows per value of ORF_name and Sequence_Name
    variants_by_ORF_by_sample = combine(groupby(df, [:ORF_name, :Sequence_Name]), nrow)

    return variants_by_ORF, variants_by_ORF_by_sample
end


variant_count427, variant_seq_count427 = count_variants(variants427)
variant_countkemp, variant_seq_countkemp = count_variants(variantskemp)