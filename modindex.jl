module ModIndex

## Compute an index modulo some given dimensions.
export modindex
modindex{N}(dims::NTuple{N, Int}, ixs::NTuple{N, Int}) = modindex(dims, ixs...)
modindex{T, N}(a::AbstractArray{T, N}, ixs::Vararg{Int, N}) = modindex(size(a), ixs...)
modindex{T, N}(a::AbstractArray{T, N}, ixs::NTuple{N, Int}) = modindex(size(a), ixs...)
@generated function modindex{N}(dims::NTuple{N, Int}, ixs::Vararg{Int, N})
    js(i) = :(mod(ixs[$i], dims[$i]))

    quote
        tuple($((:($(js(i)) == 0 ? dims[$i] : $(js(i))) for i = 1:N)...))
    end
end

## Access an element of an array at a location modulo its dimensions.
export modgetindex
modgetindex{T, N}(a::AbstractArray{T, N}, ixs::NTuple{N, Int}) = modgetindex(a, ixs...)
@generated function modgetindex{T, N}(a::AbstractArray{T, N}, ixs::Vararg{Int, N})
    js(i) = :(mod(ixs[$i], size(a, $i)))

    quote
        a[$((:($(js(i)) == 0 ? size(a, $i) : $(js(i))) for i = 1:N)...)]
    end
end

end
