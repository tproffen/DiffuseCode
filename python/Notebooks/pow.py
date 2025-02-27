import sys
#sys.path.append('/home/discus/DiffuseBuilt/python')
sys.path.append('/usr/local/lib')  # Find an installation location (/usr/local/lib or python default path)

from suite_python import suite as s
import numpy as np

s.initialize_suite()
s.execute_macro("@powder.mac")

n=200
x=np.empty(n, dtype=np.float32)
y=np.empty(n, dtype=np.float32)
s.get_data(x,y,n)

rnumbers=np.empty(4, dtype=np.float32)
s.get_r(rnumbers,1,4)

print(rnumbers)
print(x)
print(y)

s.execute_macro("@test.mac")

