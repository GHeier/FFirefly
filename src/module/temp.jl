using Firefly 

vec = Firefly.Vec(1.0, 2.0, 3.0)

function my_func(k)
    return k.x^2 + k.y^2 + k.z^2
end

println(my_func(vec))

s_val = 1.0

surf = Firefly.Surface(my_func, s_val)
@show surf.handle
println(surf)
println(surf.handle == C_NULL)
faces = get_faces(surf)
for i in 1:length(faces)
    println(faces[i][1], " ", faces[i][2])
end
