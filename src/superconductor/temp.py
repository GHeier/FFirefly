import numpy as np

wc = 0.1
T = 0.01

n = np.arange(1e7)
w_n = np.pi * T * (2 * n + 1)

term = np.where(w_n <= wc, 1 / w_n, 0)
sum = np.sum(term) * np.pi * T * 2
print(sum)
sum2 = np.sum(1 / w_n - 1 / (w_n + wc)) * 2 * np.pi * T
print(sum2)
