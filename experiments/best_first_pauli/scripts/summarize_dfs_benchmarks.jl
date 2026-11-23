python=get(ENV,"PYTHON","python")
run(`$python $(joinpath(@__DIR__,"plot_dfs_comparison.py"))`)
println("Saved DFS/BFS repeat CSV, summaries, seven comparison figures, and research report")
