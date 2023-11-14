# #input paths to folders with Variant_list.out files
folder_list = ["data/CAMP000427", "data/CAMP001339", "data/CAMP001490",
"data/CAMP001523", "data/CAMP002274", "data/CAMP003468",
"data/CAMP004884", "data/CAMP007136", "data/Kemp"] 

# load dependencies
using DataFrames
using CSV
using FileIO

#read .out files to .csv
function variant_readin(folder::AbstractString)
    # Construct the file path for Variant_list.out
    variant_file_path = joinpath(folder, "Variant_list.out")

    # Check if the file exists
    if isfile(variant_file_path)
        # Attempt to read the data from Variant_list.out with space delimiter
        try
            data = readdlm(variant_file_path, ' ', header=false)

            # Check if the file is empty or has no rows
            if isempty(data) || all(isempty, eachrow(data))
                # If empty, create an empty DataFrame with headers
                df = DataFrame(Position = String[], Original_Base = String[], Variant_Base = String[])
            else
                # Create a DataFrame
                df = DataFrame(Position = data[:, 1], Original_Base = data[:, 2], Variant_Base = data[:, 3])
            end

            # Save the DataFrame as a CSV file
            csv_file_path = joinpath(folder, "Variant_list.csv")
            CSV.write(csv_file_path, df)

            println("DataFrame saved as CSV file: $csv_file_path")

        catch
            # If reading the file fails, create an empty DataFrame with headers
            df = DataFrame(Position = String[], Original_Base = String[], Variant_Base = String[])

            # Save the empty DataFrame as a CSV file
            csv_file_path = joinpath(folder, "Variant_list.csv")
            CSV.write(csv_file_path, df, writeheader=true)

            println("Error: Unable to read Variant_list.out. Empty CSV file with headers created: $csv_file_path")
        end
    else
        println("Error: Variant_list.out not found in the specified folder.")
    end
end


for folder in folder_list
    variant_readin(folder)
end
