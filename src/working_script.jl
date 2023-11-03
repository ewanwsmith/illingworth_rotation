using FastaIO
using DataFrames
using CSV

folder_path = ["path_to_folder_1", "Path_to_folder_2"] #etc

function find_match_orfs(folder::String)
    fasta_file = joinpath(folder, "sequences.fa")
    df = fasta_readin(fasta_file)

    orfs_result = find_orfs(df)
    match_orfs_result = match_orfs(orfs_result)
    split_orfs_result = split_orfs_df(match_orfs_result)

    csv_path = joinpath(folder, "orfs.csv")
    CSV.write(csv_path, split_orfs_result)

    return split_orfs_result
end
