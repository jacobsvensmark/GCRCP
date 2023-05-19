__all__ = ["centergrid"]

import numpy as np

def centergrid(array,**kwargs):
    if ('spacing' in kwargs) and (kwargs['spacing']=='linear') or not ('spacing' in kwargs):
        array_new = (array[:-1]+array[1:])/2.
        array_new = np.insert(array_new,0,array[0]-(array_new[0]-array[0]))
        array_new = np.append(array_new,array[-1] + (array[-1]-array_new[-1]))
    elif kwargs['spacing'] == 'log':
        print("WARNING: KEYWORD LOG IS STILL JUST LINEAR")
        array_new = (array[:-1]+array[1:])/2.
        array_new = np.insert(array_new,0,array[0]-(array_new[0]-array[0]))
        array_new = np.append(array_new,array[-1] + (array[-1]-array_new[-1]))
    return array_new

