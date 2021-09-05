import xspec
from preflow import *
from lmodel_preflow import *
xspec.AllModels.addPyMod(preflow, lmodel_preflow_info, 'add')
xspec.AllData.dummyrsp(0.1, 50., 1000, scaleType='lin')
xspec.Model('preflow')
xspec.Plot.device='/xw'
xspec.Plot('emo')
