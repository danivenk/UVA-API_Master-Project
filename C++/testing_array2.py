from lagmodel.internal_class import *
import numpy as np

size = 30
shape = (6, 5)

assert size == np.prod(shape)

a = np.array(range(0, size), dtype=np.double).reshape(*shape)

print(np.array(Array(a).get_matrix()).reshape(*shape))

b = Array(*shape)

i1 = 0
i2 = 0

for i in range(0, size):
    b.set_element(i1, i2, i)
    i2 += 1

    if i1 % shape[0] == 0:
        i1 = 0
    if i2 % shape[1] == 0:
        i2 = 0
        i1 += 1

print(np.array(b.get_matrix()).reshape(*shape))
print(np.array(b))