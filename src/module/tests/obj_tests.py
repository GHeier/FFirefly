import firefly 

def point(i):
    return -1.0 + i

def val(x, y):
    return x + y

n = 3
real_file = "r_field.dat"
with open(real_file, "w") as file:
    file.write("# x y f\n")
    for i in range(n):
        for j in range(n):
            x, y = point(i), point(j)
            f = val(x, y)
            file.write(f"{x} {y} {f}\n")

real_field = firefly.Field_R(real_file)
k = [0.0, -1.0]
f = real_field(k)
print(f)


complex_file = "c_field.dat"
with open(complex_file, "w") as file:
    file.write("# x y Re(f) Im(f)\n")
    for i in range(n):
        for j in range(n):
            x, y = point(i), point(j)
            f = val(x, y)
            file.write(f"{x} {y} {f} {f}\n")

complex_field = firefly.Field_C(complex_file)
k = [0.0, -1.0]
f = complex_field(k)
print(f)
