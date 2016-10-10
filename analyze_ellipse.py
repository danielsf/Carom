import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def make_histogram(xx, dmag):
    i_xx = np.round(xx/dmag).astype(int)
    unique_ixx, ct = np.unique(i_xx, return_counts=True)

    return unique_ixx*dmag, ct.astype(float)


data_file = 'ellipse_samples.txt'
control_file = 'junk.txt'
control_dtype=np.dtype([('dim', int), ('step', int), ('pdf', float)])
control_data = np.genfromtxt(control_file, control_dtype)

dtype = np.dtype([('x1', float), ('x2', float), ('x3', float), ('x4', float)])

radii= {}
radii['x1']=0.3539757
radii['x2']=1.438721
radii['x3']=1.724201
radii['x4']=1.501453

data = np.genfromtxt(data_file, dtype=dtype)

dx=0.1

plt.figsize=(30,30)
plt.subplot(1,2,1)


plt.plot(control_data['step'].astype(float)/float(control_data['step'].max()),
         control_data['pdf'], color='k')

for col, color in zip(('x1', 'x2', 'x3', 'x4'),
                      ('b', 'r', 'g', 'm')):
    
    samples = np.abs(data[col]/radii[col])
    xvals, ct = make_histogram(samples, dx)
    ct = ct/ct.max()
    ct = ct*control_data['pdf'][0]
    plt.plot(xvals, ct,color=color)


plt.subplot(1,2,2)
dtype = np.dtype([('name', str, 10), ('roll', float), ('other_name', str, 10),
                  ('dex', int), ('third_name', str, 10), ('coord', float)])

check_data = np.genfromtxt('other_junk.txt', dtype=dtype)
samples=np.abs(check_data['dex'])
vals, ct = make_histogram(samples, 0.1)
plt.plot(vals, ct,color='k')

plt.savefig('ellipse_plot.png')
    
    
    
