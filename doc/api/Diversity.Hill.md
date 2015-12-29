# Diversity.Hill


## Methods [Exported]

---

<a id="method__hillnumber.1" class="lexicon_definition"></a>
#### hillnumber(proportions,  qs) [¶](#method__hillnumber.1)
### Calculates Hill numbers

Calculate the Hill number (or naive diversity) of order q of
population(s) with given relative proportions

### Arguments:
- `proportions`: relative proportions of different individuals / species
                 in population (vector, or matrix where columns are
                 individual populations) 
- `qs`: single number or vector of orders of diversity measurement

### Returns:
- Diversity of order qs (single number or vector of diversities)


*source:*
[Diversity/src/Hill.jl:18](https://github.com/richardreeve/Diversity.jl/tree/33e5ccdee2b395af7d54b11ad13cb4d67c8531ce/src/Hill.jl#L18)

