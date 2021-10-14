# function to write logarithmic x-data to a file
a1 = np.array([1,2,3,4,5,6,7,8,9])
a2 = np.array([1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9])
a3 = np.sort(np.outer(a2, a1).flatten())

a=""
for i in a3:
    a=a+str(i)+", "

f = open("x_data.txt", "w")
f.write(a)
f.close()