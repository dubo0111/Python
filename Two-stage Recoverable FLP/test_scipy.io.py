import scipy.io
mat = scipy.io.loadmat('city_250.mat')
a=mat['city_250']
#print(a)
b=a[0,:]
#print(b)
scipy.io.savemat('save_b.mat',{'b':b})
bb = scipy.io.loadmat('save_b.mat')

