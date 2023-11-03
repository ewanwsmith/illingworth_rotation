using DataFrames
using CSV

#input paths to folders with .fasta files
folder_list = ["data/CAMP000427", "data/CAMP001339", "data/CAMP001490",
                "data/CAMP001523", "data/CAMP002274", "data/CAMP003468",
                "data/CAMP004884", "data/CAMP007136", "data/Kemp"] 
#input path to reference ORFs
reference_path = "data/reference/coding_sequences.fasta" 

function find_match_orfs(folder::String)
    ref_orfs = ref_readin(reference_path)
    fasta_file = joinpath(folder, "sequences.fa")
    df = fasta_readin(fasta_file)

    orfs_result = find_orfs(df)
    match_orfs_result = match_orfs(orfs_result, ref_orfs)

    csv_path = joinpath(folder, "matched_orfs.csv")
    CSV.write(csv_path, match_orfs_result)

    return match_orfs_result
end

for folder in folder_list
    find_match_orfs(folder)
end