# SPDX-License-Identifier: BSD-2-Clause

module DiversityAxisArraysExt

isdefined(Base, :get_extension) ? (using AxisArrays) : (using ..AxisArrays)
import Diversity

function Diversity.GeneralTypes(zmatrix::AM) where
         {FP <: AbstractFloat, M <: AbstractMatrix{FP}, LR, LC, NAMES,
          T <: Tuple{Axis{LR, NAMES}, Axis{LC, NAMES}},
          AM <: AxisMatrix{FP, M, T}}
    size(zmatrix, 1) == size(zmatrix, 2) ||
        throw(DimensionMismatch("Similarity matrix is not square"))

    AxisArrays.axes(zmatrix, 1).val == AxisArrays.axes(zmatrix, 2).val ||
        throw(DimensionMismatch("Similarity matrix does not have matching row and column labels"))

    return Diversity.GeneralTypes(zmatrix,
                                  collect(AxisArrays.axes(zmatrix, 1).val))
end

end
