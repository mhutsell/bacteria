#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 17:06:43 2020

@author: Mack
"""

def diff_normal(y):
    n, gmax, row, k, mu, a = y
    dndt = (n*gmax*row)/(row+k) - (n*mu)
    dpdt = (-a)*(n*gmax*row)/(row+k)
    diff = [dndt, dpdt]

def diff_pre(y):
    n, gmax, row, k, mu, a = y
    dndt = -(n*mu)
    dpdt = 0
    diff = [dndt, dpdt]