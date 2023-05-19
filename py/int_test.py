import numpy as np
from scipy import interpolate

x = np.array([1, 2, 3])
y = np.array([4, 5])
xx, yy = np.meshgrid(x, y)

z = np.array(
    [
        [0.0, 5.0, 4.0],
        [0.0, 10.0, 8.0],
    ]
)

print("x")
print(x)
print("y")
print(y)
print("xx")
print(xx)
print("yy")
print(yy)
print("z")
print(z)
print("np.shape(z)")
print(np.shape(z))

print(np.shape(x))
print(np.shape(y))
print(np.shape(z))

f_1 = interpolate.interp2d(xx, yy, z, kind="linear")
f_2 = interpolate.interp2d(x, y, z, kind="linear")
f_3 = interpolate.interp2d(np.ravel(xx), np.ravel(yy), np.ravel(z), kind="linear")

xnew = np.array([1, 2.1, 3])
ynew = np.array([4, 5, 6])

znew_1 = f_1(xnew, ynew)
znew_2 = f_2(xnew, ynew)
znew_3 = f_3(xnew, ynew)

