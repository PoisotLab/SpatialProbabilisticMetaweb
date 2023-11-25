files = readdir("./data/results/motifs/")

S1 = filter(contains("S1"), files)
ids = parse.(Int, replace.(S1, "S1-" => "", ".tif" => ""))
[setdiff(1:500, ids)]

S2 = filter(contains("S2"), files)
ids = parse.(Int, replace.(S2, "S2-" => "", ".tif" => ""))
[setdiff(1:500, ids)]

S4 = filter(contains("S4"), files)
ids = parse.(Int, replace.(S4, "S4-" => "", ".tif" => ""))
setdiff(1:500, ids) # 245

S5 = filter(contains("S5"), files)
ids = parse.(Int, replace.(S5, "S5-" => "", ".tif" => ""))
setdiff(1:500, ids) # 51, 86, 270, 427, 431, 463

files = readdir("jobs/out/"; join=true)
filter!(contains("13_get_motifs_S4"), files)
