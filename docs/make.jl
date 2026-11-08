using Documenter
using Literate
using Plots
using GravitationalWaveData

# Show Plots.jl figures in the built docs as PNG files. Documenter uses the first of text/html,
# SVG and PNG that an @example result offers; a Plots.jl figure offers all three, and for strain
# plots with tens of thousands of samples the HTML/SVG forms run to megabytes each. Declining the
# first two for GR figures here makes Documenter write a compact PNG instead. This only affects
# the docs build: the tutorials and notebooks themselves are unchanged.
Base.showable(::MIME"text/html", ::Plots.Plot{Plots.GRBackend}) = false
Base.showable(::MIME"image/svg+xml", ::Plots.Plot{Plots.GRBackend}) = false

const DOCS_SRC = joinpath(@__DIR__, "src")
const LITERATE_DIR = joinpath(@__DIR__, "literate")
const GENERATED_DIR = joinpath(DOCS_SRC, "generated")

# Weave every docs/literate/*.jl tutorial into docs/src/generated/*.md (build
# artifact, gitignored). The same sources produce notebooks/*.ipynb via
# notebooks/generate.jl.
# Ordered explicitly (rather than a directory listing) so the sidebar order in
# `pages` below matches; add new topics to BOTH this list and `pages`.
const LITERATE_TOPICS = [
    "01_getting_started",
    "02_psd_asd",
    "03_filtering_whitening",
    "04_resampling",
    "05_spectrogram",
    "06_qtransform",
    "07_gps_segments",
    "08_detector_response",
    "09_coherence",
    "10_gwosc_gwf_fetch",
    "11_plotting_recipes",
    "12_gw150914_capstone",
]

# docs/CONVENTIONS.md (repo root) is the single source of truth; copy it in
# verbatim rather than maintaining a second copy under docs/src/.
cp(joinpath(@__DIR__, "CONVENTIONS.md"), joinpath(DOCS_SRC, "conventions.md"); force=true)

mkpath(GENERATED_DIR)
for topic in LITERATE_TOPICS
    Literate.markdown(joinpath(LITERATE_DIR, "$(topic).jl"), GENERATED_DIR; documenter=true)
end

generated_pages = [
    "Getting Started" => "generated/01_getting_started.md",
    "PSD & ASD Estimation" => "generated/02_psd_asd.md",
    "Filtering & Whitening" => "generated/03_filtering_whitening.md",
    "Resampling" => "generated/04_resampling.md",
    "Spectrogram" => "generated/05_spectrogram.md",
    "Q-transform" => "generated/06_qtransform.md",
    "GPS Time & Segments" => "generated/07_gps_segments.md",
    "Detector Response" => "generated/08_detector_response.md",
    "Coherence & CSD" => "generated/09_coherence.md",
    "GWOSC & GWF Fetch" => "generated/10_gwosc_gwf_fetch.md",
    "Plotting Recipes" => "generated/11_plotting_recipes.md",
    "Capstone: GW150914" => "generated/12_gw150914_capstone.md",
]

DocMeta.setdocmeta!(GravitationalWaveData, :DocTestSetup, :(using GravitationalWaveData); recursive=true)

makedocs(;
    sitename="GravitationalWaveData.jl",
    modules=[GravitationalWaveData],
    authors="Paul Hopkins <paul.hopkins256@gmail.com>",
    format=Documenter.HTML(;
        canonical="https://juliagravitationalwaves.github.io/GravitationalWaveData.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => generated_pages,
        "Conventions" => "conventions.md",
        "API Reference" => [
            "Overview" => "api/index.md",
            "Data types & I/O" => "api/data.md",
            "Time & Segments" => "api/time.md",
            "Spectral Analysis" => "api/spectral.md",
            "Filtering & Resampling" => "api/processing.md",
            "Detectors & Reference Curves" => "api/detectors.md",
        ],
    ],
    checkdocs=:exports,
    warnonly=[:missing_docs, :cross_references],
)

deploydocs(;
    repo="github.com/JuliaGravitationalWaves/GravitationalWaveData.jl",
    devbranch="main",
)
