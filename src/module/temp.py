import firefly
import matplotlib.pyplot as plt

vec = firefly.Vec(1.0, 2.0, 3.0)


def my_func(k):
    return k.x**2 + k.y**2 + k.z**2


s_val = 1.0

surf = firefly.Surface(my_func, s_val)
faces = surf.faces
print(faces[0], faces[0][0])
print(faces[1], faces[1][0])
print(faces[2], faces[2][0])
print(faces[3], faces[3][0])

plt.figure()
for i in range(len(faces)):
    plt.scatter(faces[i][0], faces[i][1])
plt.show()
