import sys
sys.path.append('/home/discus/DiffuseBuilt/python')
#sys.path.append('/usr/local/lib')  # Find an installation location (/usr/local/lib or python default path)

from suite_python import suite as s
import numpy as np

s.initialize_suite()
s.execute_macro("@one.mac")
s.execute_macro("@two.mac")
