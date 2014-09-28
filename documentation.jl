using Lexicon
using Diversity
using Diversity.Ecology
using Diversity.Hill
using Diversity.Jost

function createdocs(dir::String)
    slash = (dir[end] == '/') ? "" : "/"
    save("$(dir)$(slash)diversity.html", Diversity)
    save("$(dir)$(slash)ecology.html", Diversity.Ecology)
    save("$(dir)$(slash)hill.html", Diversity.Hill)
    save("$(dir)$(slash)jost.html", Diversity.Jost)
end
