# Get parameters from excel

using XLSX, DataFrames

# open .xlsx file
xlsx_file = XLSX.readxlsx("./jl_code/exampleNet.xlsx")
# print sheetnames
sheet_name = XLSX.sheetnames(xlsx_file)

# species, reactions
species_sheet = xlsx_file[sheet_name[1]]
reactions_sheet = xlsx_file[sheet_name[2]]

# change to DataFrame
df_species_raw = DataFrame(species_sheet[:], :auto)
df_reactions_raw = DataFrame(reactions_sheet[:], :auto)

# get # of species
column_ID = df_species_raw[:, 2]
ID = filter(!ismissing, column_ID)
num_species = size(ID[2:end, :]) 

y0 = df_species_raw[2:num_species]
tau = xlsx_file["C:C"]

y0 = hvcat(size(y0,2), permutedims(y0)...)
ymax = hvcat(size(ymax,2), permutedims(ymax)...)
tau = hvcat(size(tau,2), permutedims(tau)...)

w = xlsx_file["reactions!D3:D8"]
n = xlsx_file["reactions!E3:E7"]
EC50 = xlsx_file["reactions!F3:F7"]
        
w = hvcat(size(w,2), permutedims(w)...)
n = hvcat(size(n,2), permutedims(n)...)
EC50 = hvcat(size(EC50 ,2), permutedims(EC50)...)