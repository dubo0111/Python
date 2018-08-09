import numpy as np
import matplotlib.pyplot as plt
from plotly import offline as py
py.init_notebook_mode()
t = np.linspace(0, 20, 500)
plt.plot(t, np.sin(t))
py.iplot_mpl(plt.gcf())
from IPython.display import Math
Math(r'i\hbar \frac{dA}{dt}~=~[A(t),H(t)]+i\hbar \frac{\partial A}{\partial t}.')
import pandas as pd
df = pd.DataFrame({'A': 1.,
                   'B': pd.Timestamp('20130102'),
                   'C': pd.Series(1, index=list(range(4)), dtype='float32'),
                   'D': np.array([3] * 4, dtype='int32'),
                   'E': pd.Categorical(["test", "train", "test", "train"]),
                   'F': 'foo'})

df
