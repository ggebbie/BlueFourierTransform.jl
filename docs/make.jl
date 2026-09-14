using BlueFourierTransform
using Documenter

DocMeta.setdocmeta!(BlueFourierTransform, :DocTestSetup, :(using BlueFourierTransform); recursive=true)

makedocs(;
    modules=[BlueFourierTransform],
    authors="G Jake Gebbie <ggebbie@whoi.edu>",
    sitename="BlueFourierTransform.jl",
    format=Documenter.HTML(;
        canonical="https://ggebbie.github.io/BlueFourierTransform.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ggebbie/BlueFourierTransform.jl",
    devbranch="main",
)
