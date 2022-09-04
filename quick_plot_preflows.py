import xspec
from preflows import *
from lmodel_preflow import *
xspec.AllModels.addPyMod(preflows, lmodel_preflows_info, 'add')
xspec.AllData.dummyrsp(0.01, 50., 1000, scaleType='log')
xspec.Model('preflows')
xspec.Plot.device='/xw'
xspec.Plot('emo')
