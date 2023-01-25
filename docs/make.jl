push!(LOAD_PATH,"../src/")
using sample_project, Documenter
makedocs(
         sitename = "sample_project.jl"
         , 
        #  modules  = [HomogeneityTestBBU],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(
    repo="github.com/bunting-econ/sample_project.jl.git",
)
