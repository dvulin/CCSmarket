# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 11:58:12 2022

@author: dvulin
"""

"""
klase:
    emisije
        #beccs, carbon farming
        #industry (cement, fertilizer)
        #storageOperator
        #service
        
"""

class Emissions:
    import numpy as np
    def __init__(self):
        self.test = 'pass'
          
    def addTimeSeries(self, time_steps):
        """
        Parameters
        ----------
        time_steps : dodaje listu dana, datuma i sl., ovisno o tipu varijable

        Returns
        -------
        None.

        """
        if type(time_steps) == 'int':
            time_steps = np.arange(time_steps)
            
class Industry(Emissions):
    def __init__(self):
        pass
    
         
        

from emissions import Emissions
e = Emissions()


