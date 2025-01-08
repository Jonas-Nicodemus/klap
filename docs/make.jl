CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
using DrWatson
@quickactivate "klap"
using Documenter, DocumenterCitations
using Passivation

bib = CitationBibliography(joinpath(@__DIR__, "..", "CITATION.bib"))

@info "Building Documentation"
makedocs(;
    sitename = "klap",
    # This argument is only so that the sequence of pages in the sidebar is configured
    # By default all markdown files in `docs/src` are expanded and included.
    pages=[
        "KLAP" => "index.md",
        "Passivation" => "Passivation.md",
        "API" => "API.md",
    ],
    # Don't worry about what `CI` does in this line.
    format=Documenter.HTML(;
        prettyurls=CI,
        canonical="https://Jonas-Nicodemus.github.io/klap",
        edit_link="main",
        assets=String[],
    ),
    plugins=[bib],
)

@info "Deploying Documentation"
if CI
    deploydocs(
        # `repo` MUST be set correctly. Once your GitHub name is set
        # the auto-generated documentation will be hosted at:
        # https://Jonas-Nicodemus.github.io/klap/dev/
        # (assuming you have enabled `gh-pages` deployment)
        repo = "github.com/Jonas-Nicodemus/klap.git",
        target = "build",
        push_preview = true,
        devbranch = "main",
    )
end

@info "Finished with Documentation"
