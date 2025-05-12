using Firefly 

vec = Firefly.Vec(1.0, 2.0, 3.0)
println(vec.x)
println(vec.y)
println(vec.z)

function my_func(k)
    println("called")
    println(k.x, k.y, k.z)
    return k.x^2 + k.y^2 + k.z^2
end

println(my_func(vec))

s_val = 1.0

surf = Firefly.Surface(my_func, s_val)
@show surf.handle
println(surf)
println(surf.handle == C_NULL)
faces = get_faces(surf)
println(faces)
