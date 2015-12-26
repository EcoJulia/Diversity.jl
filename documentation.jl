using Lexicon
using Docile
using Diversity
using Diversity.Ecology
using Diversity.Hill
using Diversity.Jost
@docstrings

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
    slash = (dir[end] == '/') ? "" : "/"
    save("$(dir)$(slash)diversity.html", Diversity)
    save("$(dir)$(slash)ecology.html", Diversity.Ecology)
    save("$(dir)$(slash)hill.html", Diversity.Hill)
    save("$(dir)$(slash)jost.html", Diversity.Jost)
end

"""
### createmddocs()

Create markdown documentation for Diversity, Diversity.Ecology,
Diversity.Hill and Diversity.Jost packages in the directory given. The
documentation is currently stored in \"doc/api\". These will then be
uploaded to github.
"""
function createmddocs(dir::AbstractString)
    slash = (dir[end] == '/') ? "" : "/"
    save("$(dir)$(slash)Diversity.md", Diversity)
    save("$(dir)$(slash)Diversity.Ecology.md", Diversity.Ecology)
    save("$(dir)$(slash)Diversity.Hill.md", Diversity.Hill)
    save("$(dir)$(slash)Diversity.Jost.md", Diversity.Jost)
    index = Index()
    update!(index, save("doc/api/Diversity.md", Diversity));
    update!(index, save("doc/api/Diversity.Ecology.md", Diversity.Ecology));
    update!(index, save("doc/api/Diversity.Hill.md", Diversity.Hill));
    update!(index, save("doc/api/Diversity.Jost.md", Diversity.Jost));
    save("$(dir)$(slash)index.md", index)
end
