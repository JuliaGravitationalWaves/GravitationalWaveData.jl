# GWData.jl notebooks

- `*.ipynb` — Jupyter notebooks, one per docs tutorial, generated from `docs/literate/*.jl` via
  `julia --project=docs generate.jl`. Do not hand-edit; edit the Literate source under
  `docs/literate/` instead and regenerate. Requires an `IJulia` kernel to run (not a dependency of
  this package itself — install it in your own global environment: `using Pkg; Pkg.add("IJulia")`).
- `pluto/*.jl` — standalone interactive [Pluto.jl](https://plutojl.org/) notebooks (not generated;
  edited directly). Open with `using Pluto; Pluto.run()` then open the file.
- `bonito/*.jl` — [Bonito.jl](https://bonito.sh/) versions of the same interactive explorers as
  `pluto/*.jl` (not generated; edited directly). Each is a standalone script that starts its own
  web server rather than needing a Pluto/Jupyter frontend: run e.g.
  `julia --project=notebooks/bonito notebooks/bonito/filter_explorer.jl` and open the printed
  `http://127.0.0.1:PORT` URL in a browser (port 8081 for the filter explorer, 8082 for the
  Q-transform explorer). The first page load takes a minute or so while the data loads and the
  first plot renders. Uses its own environment (`bonito/Project.toml`, with GWData taken from the
  repo root via `[sources]`, so Julia 1.11 or later) rather than `docs`, because Bonito 5.x needs
  HTTP.jl 2.x while Pluto pins HTTP.jl to 1.x.
