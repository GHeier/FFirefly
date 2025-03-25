module FieldTests

include("../field.jl")
using .CondensedMatterField

using Test

function pretty_print(mesh)
    sz, dim = size(mesh)
    # Print points
    for i in 1:sz
        print("(", mesh[i, 1])
        for j in 2:dim
            print(", ", mesh[i, j])
        end
        println(")")
    end
end

function build_mesh(point_list)
    dim = length(point_list)
    sz = length(point_list[1])
    points = Matrix{Float64}(undef, sz^dim, dim)
    for i in 1:dim, j in 1:sz, k in 1:sz
        ind = (j - 1) * sz + k
        ind1 = Int(trunc(ind / sz)) + 1
        if ind1 > sz
            ind1 = sz
        end
        ind2 = ind % sz
        if ind2 == 0
            ind2 = sz
        end
        points[ind, i] = point_list[i][i%2 * ind1 + (i-1)%2 * ind2]
    end
    println("size: ", size(points))
    #pretty_print(points)
    return points
end

function test_2d_field()
    # Generate points
    x = range(-1, 1, length=10)
    y = range(-1, 1, length=10)

    x = Vector(x)
    y = Vector(y)
    combined = [x, y]
    mesh = build_mesh(combined)
    data = x .^ 2 .+ y .^ 2
    a = CMF(mesh, data, 2, false)

    r1 = a(0.5, 0.5)
    r2 = a(1, 1)
    r3 = a(0, 0)
    r4 = a(-1, -1)

    @test r1 ≈ 0.5^2 + 0.5^2
    @test r2 ≈ 1^2 + 1^2
    @test r3 ≈ 0^2 + 0^2
    @test r4 ≈ (-1)^2 + (-1)^2
end

end # module

using .FieldTests
FieldTests.test_2d_field()
