#!/usr/bin/env python3
from skl2onnx import to_onnx
import  numpy as np
import pandas as pd
import pickle
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


modelDir = "models"
sels = [
        "CR+2J",
        "WALP",
        "WHS_highET",
        "WHS_lowET",
        "ZHS_highET",
        "ZHS_lowET",
        ]

datasets = ["ggH125_S55_ctau5p32_gg_py8.csv", "mALP1_W_py8_ct0.031.csv", "WH400_S100_W_py8_ct1.25.csv" , "WH125_S16_W_py8_ct0.3.csv", "mHZZd_600_150_Z_py8_ct1.6.csv", "mHZZd_250_50_Z_py8_ct3.4.csv"]

for i, sel in enumerate(sels):
   df = pd.read_csv(datasets[i])
   with open(f"{modelDir}/{sel}_features.txt") as g: var=eval(g.read())
   f = open(f'{modelDir}/{sel}_model.pkl', 'rb')
   mean =  np.load(f"{modelDir}/{sel}_scaler_mean.npy")
   std =  np.load(f"{modelDir}/{sel}_scaler_std.npy")
   print(mean)
   print(std)
   #var = [v.replace("W","V").replace("Z","V") for v in var] 
   scaler = StandardScaler() 
   scaler.mean_ = mean
   scaler.scale_ = np.ones(len(mean))
   scaler.scale_ = std
   clf = pickle.load(f)
   model = Pipeline(steps=[
   ('scaler', scaler),
   ('classifier', clf)
   ])
   #inData= [[2.48167911, 1.54975031, 0.58969972, 156.00026434, 156.32044804, 21., 2.48167911, 1.54975031, 0.58969972, 156.00026434, 156.32044804, 21. , 156.8854222 , -0.98841971]]
   #print (list(inData))
   #print(list(scaler.transform(inData)))
   #exit(1)
   #X = model.predict_proba(inData)
   #print(X)
   #onx = to_onnx(clf, {var[i]: df[var].iloc[0].to_numpy()[i] for i in range(len(var))})
   inputsize = len(var)
   from skl2onnx.common.data_types import FloatTensorType, Int64TensorType
   final_type = [
               ('output', FloatTensorType([1, 1])),
               ('output_probability', FloatTensorType([1, 5])),
                 ]
   onx = to_onnx(model,  initial_types=[('float_input', FloatTensorType([1, inputsize]))] , final_types=final_type,  options={'zipmap': False})
   #onx = to_onnx(model,  df[var].values[:1], initial_types=[('float_input', FloatTensorType([1, inputsize]))])
   with open(f"{modelDir}/{sel}_model.onnx", "wb") as h:
     h.write(onx.SerializeToString())

   #import onnx.shape_inference

   ## Load and infer shapes
   #model = onnx.load(f"{modelDir}/{sel}_model.onnx")
   #inferred_model = onnx.shape_inference.infer_shapes(model)

   ## Print inferred types for all outputs
   #for output in inferred_model.graph.output:
   #  print(f"Output Name: {output.name}")
   #  print(f"Output Type: {output.type}")

