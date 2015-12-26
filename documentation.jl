using Lexicon, Docile
using Diversity, Diversity.Ecology, Diversity.Hill, Diversity.Jost

"""
### createhtmldocs()

Create html documentation for Diversity, Diversity.Ecology,
Diversity.Hill and Diversity.Jost packages in the directory given. The
master documentation is currently stored in \"doc/site/master\", and
the stable documentation for the latest released version of the
packages is stored in \"doc/site/stable\". These will then be uploaded
to github and displayed on github.io using:

git subtree push --prefix doc/site origin gh-pages
"""
function createhtmldocs(dir::AbstractString)
    save(joinpath(dir, "diversity.html"), Diversity)
    save(joinpath(dir, "ecology.html"), Diversity.Ecology)
    save(joinpath(dir, "hill.html"), Diversity.Hill)
    save(joinpath(dir, "jost.html"), Diversity.Jost)
end

"""
### createmddocs()

Create markdown documentation for Diversity, Diversity.Ecology,
Diversity.Hill and Diversity.Jost packages in the directory given. The
documentation is currently stored in \"doc/api\". These will then be
uploaded to github.
"""
function createmddocs(dir::AbstractString)
    const modules = Docile.Collector.submodules(Diversity)
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
        update!(index, save(joinpath(dir, "$(mod).md"), mod, config))
    end
    save(joinpath(dir, "index.md"), index, config)
    
    # Add a reminder not to edit the generated files.
    open(joinpath(dir, "README.md"), "w") do f
        print(f, """
        Files in this directory are generated using an automated script. Make
        all changes to the originating docstrings/files rather than these ones.

        Documentation should *only* be built directly on the `master` branch.
        Source links would otherwise become unavailable should a branch be
        deleted from the `origin`. This means potential pull request authors
        *should not* run the build script when filing a PR.
        """)
    end

    info("Need to build website using mkdocs from Diversity root:")
    info("Run \"mkdocs build --clean\"")
end
