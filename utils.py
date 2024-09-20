#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 12:46:30 2023

@author: morgan.schneider

Helpful all-purpose functions and classes
"""

####################
### Load modules ###
####################




#%%
######################
### Define classes ###
######################

class struct():
    def __init__(self,**kwargs):
        self.Set(**kwargs)
    def Set(self,**kwargs):
        self.__dict__.update(**kwargs)
    def SetAttr(self,var,val):
        self.__dict__[var] = val
    def keys(self):
        print(self.__dict__.keys())