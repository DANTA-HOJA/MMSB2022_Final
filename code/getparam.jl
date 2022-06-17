# Get parameters from excel
using XLSX

#y0 = XLSX.readdata("exampleNet.xlsx", "species", "D3:D7")
raw = XLSX.readxlsx("exampleNet.xlsx")
y0 = raw["species!D3:D7"]
ymax = raw["species!E3:E7"]
tau = raw["species!F3:F7"]

y0 = hvcat(size(y0,2), permutedims(y0)...)
ymax = hvcat(size(ymax,2), permutedims(ymax)...)
tau = hvcat(size(tau,2), permutedims(tau)...)

w = raw["reactions!D3:D8"]
n = raw["reactions!E3:E7"]
EC50 = raw["reactions!F3:F7"]
        
w = hvcat(size(w,2), permutedims(w)...)
n = hvcat(size(n,2), permutedims(n)...)
EC50 = hvcat(size(EC50 ,2), permutedims(EC50)...)

reactionIDs = raw["reactions!B3:B8"]
reactionRules = raw["reactions!C3:C8"]
specID = raw["species!B2:B7"]


