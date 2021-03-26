# -*- coding: utf-8 -*-
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2020/04/10*

*licence : AGPL-3.0*

Module to get screen resolution 
    ..
"""

import tkinter as tk

root_tmp = tk.Tk()
screen_width = root_tmp.winfo_screenwidth()
screen_height = root_tmp.winfo_screenheight()
print("Screen resolution")
print('width       :',screen_width)
print('height      :',screen_height)
