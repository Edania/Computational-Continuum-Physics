import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.insert(0, 'C:/Users/elsa2/Documents/Master/Computational Continuum Physics/test/build/Release')
import testing

#TODO: Implement residual, test for convergence, approximation order and so on
print('nyeh')
result = testing.advection(100,100, 0.001, 0.001)

x = np.linspace(0,100,100).reshape(1,-1)
t = np.linspace(0,100,100).reshape(-1,1)
analytical = x-t

print(result.shape)
print(analytical.shape)

fig = plt.figure()
plt.imshow(result)
plt.title('numerical')
plt.xlabel('x')
plt.ylabel('t')
fig.savefig('numerical.pdf')

fig = plt.figure()
plt.imshow(x-t)
plt.title('analytical')
plt.xlabel('x')
plt.ylabel('t')
fig.savefig('analytical.pdf')

