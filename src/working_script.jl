using FastaIO
using DataFrames
using CSV

folder_list = ["path_to_folder_1", "Path_to_folder_2"] #input paths to folders with .fasta files
reference_path = "path_to_reference_ORFs" #input paths to reference ORFs

function find_match_orfs(folder::String)
    fasta_file = joinpath(folder, "sequences.fa")
    df = fasta_readin(fasta_file)

    orfs_result = find_orfs(df)
    match_orfs_result = match_orfs(orfs_result)

    csv_path = joinpath(folder, "matched_orfs.csv")
    CSV.write(csv_path, match_orfs_result)

    return match_orfs_result
end

for folder in folder_list
    find_match_orfs(folder)
end