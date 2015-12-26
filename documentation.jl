using Lexicon
using Docile
using Diversity
using Diversity.Ecology
using Diversity.Hill
using Diversity.Jost
@docstrings

@doc """
### createdocs()

Create html documentation for Diversity, Diversity.Ecology,
Diversity.Hill and Diversity.Jost packages in the directory given. The
master documentation is currently stored in \"doc/site/master\", and
the stable documentation for the latest released version of the
packages is stored in \"doc/site/stable\". These will then be uploaded
to github and displayed on github.io using:

git subtree push --prefix doc/site origin gh-pages
""" ->
function createdocs(dir::AbstractString)
    slash = (dir[end] == '/') ? "" : "/"
    save("$(dir)$(slash)diversity.html", Diversity)
    save("$(dir)$(slash)ecology.html", Diversity.Ecology)
    save("$(dir)$(slash)hill.html", Diversity.Hill)
    save("$(dir)$(slash)jost.html", Diversity.Jost)
end
