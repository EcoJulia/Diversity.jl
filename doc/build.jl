using Diversity, Diversity.Ecology, Diversity.Hill, Diversity.Jost
using Docile, Lexicon

const api_directory = "api"
const modules = Docile.Collector.submodules(Diversity)

cd(dirname(@__FILE__)) do

    # Run the doctests *before* we start to generate *any* documentation.
    for m in modules
        failures = failed(doctest(m))
        if !isempty(failures.results)
            println("\nDoctests failed, aborting commit.\n")
            display(failures)
            exit(1) # Bail when doctests fail.
        end
    end

    # Generate and save the contents of docstrings as markdown files.
    index  = Index()
    config = Config(md_subheader = :category,
                    category_order = [:module, :function, :method,
                                      :type, :typealias,
                                      :macro, :global])
    for mod in modules
        update!(index, save(joinpath(api_directory, "$(mod).md"), mod, config))
    end
    save(joinpath(api_directory, "index.md"), index, config)

    # Add a reminder not to edit the generated files.
    open(joinpath(api_directory, "README.md"), "w") do f
        print(f, """
        Files in this directory are generated using the `build.jl` script. Make
        all changes to the originating docstrings/files rather than these ones.

        Documentation should *only* be build directly on the `master` branch.
        Source links would otherwise become unavailable should a branch be
        deleted from the `origin`. This means potential pull request authors
        *should not* run the build script when filing a PR.
        """)
    end

    info("Need to run \"git add $(api_directory)\" to add new api files from this directory.")
end
