# add equations and meta-code to methods

import math
import numpy as np
from scipy.stats import gmean, hmean

p = [0.000000021, 0.00417922, 0.000131232, 0.17]
p1 = [0.021, 0.00417922, 0.000131232, 0.17]
p2 = [0.21, 0.417922, 0.0131232, 0.17]
p3 = [0.0021, 0.00417922, 0.00131232, 0.17]

math.e * gmean(p)
math.e * gmean(p1)
math.e * gmean(p2)
math.e * gmean(p3)

math.e * np.log(len(p)) * hmean(p)
math.e * np.log(len(p1)) * hmean(p1)
math.e * np.log(len(p2)) * hmean(p2)
math.e * np.log(len(p3)) * hmean(p3)





def median(lst):
    n = len(lst)
    s = sorted(lst)
    return (sum(s[n//2-1:n//2+1])/2.0, s[n//2])[n % 2] if n else None

def mean(lst):
    return (sum(lst) / len(lst))


l = [-2.128803615, -2.961825969, -0.494695342, 0.38]
l1 = [-2.128803615, -2.961825969, -2.494695342, 0.38]
l2 = [-2.128803615, -2.961825969, 0.494695342, 0.38]
l3 = [-1.128803615, -1.961825969, -0.494695342, 0.38]

median(l)
median(l1)
median(l2)
median(l3)

mean(l)
mean(l1)
mean(l2)
mean(l3)
