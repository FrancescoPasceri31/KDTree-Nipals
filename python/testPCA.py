import numpy as np

np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)})

vet = [[-28.0, -75.0, -140.0, -52.0], [-14.0, -33.0, -100.0, -36.0],
       [-9.0, -21.0, -97.0, -27.0], [-3.0, 2.0, -1.0, -15.0],
       [1.0, 5.0, 10.0, -9.0], [5.0, 23.0, 11.0, -2.0],
       [27.0, 61.0, 15.0, 12.0], [32.0, 172.0, 27.0, 24.0],
       [33.0, 340.0, 40.0, 36.0]]
n = 9
k = 4
ds = np.array(vet)
#print(ds)
#print(np.transpose(ds))

teta = 1 * np.exp(-8)
print()
#print(teta)

for i in range(0, k):
    mean = 0.0
    for j in range(0, n):
        mean += ds[j][i]
    mean /= n
    #print(mean)
    for j in range(0, n):
        ds[j][i] -= mean
#print()
#print(ds)
#print(np.shape(ds))
#print()

h = 2
U = np.zeros(((n),(h)))
V = np.zeros(((k),(h)))
u = np.transpose(np.array([ ds[:, 0] ]))


for j in range(0, h):
    while True:

        x = (u.T).dot(u)[0][0]
        v = np.transpose(ds).dot(u) / x
        v = v / np.linalg.norm(v)

        t = x

        u = ds.dot(v) / (v.T).dot(v)
        
        tt = (u.T).dot(u)[0][0]

        if(np.abs(tt - t) < teta * tt):
            break
    U[:,j] = u[:,0]
    V[:,j] = v[:,0]
    ds = ds - np.dot(u,np.transpose(v))
print()
print(U)
print()
print(V)
print()
print(U.dot(V.T))
