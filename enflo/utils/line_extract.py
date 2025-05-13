#!/usr/bin/env python
"""
Usage: 
Run "line_extract x0 x1 Nx <vtk-file>" to extract the line data from particular vtk file
Run "line_extract x0 x1 Nx" to extract the line data all vtk files of the form sol**.vtk

"""

import sys,glob
import math
import vtk
import matplotlib.pyplot as plt
import numpy as np

#Check if a vtk file argument is given
if len(sys.argv) < 4:
   print("Please specify x0,x1,Nx, vtkfile(optional).")
   print("Example:")
   print(f"  line_extract.py 0 1 256 sol200.vtk")
   print("If vtk file is not specified, all files with the name sol*.vtk will be read. ")
   sys.exit(1)
elif len(sys.argv)==4:
   x0 = float(sys.argv[1])
   x1 = float(sys.argv[2])
   Nx = int(sys.argv[3])
   vtkfiles = glob.glob("sol*.vtk")
   print(vtkfiles)

elif len(sys.argv)==5:
   x0 = float(sys.argv[1])
   x1 = float(sys.argv[2])
   Nx = int(sys.argv[3])
   vtkfiles = [sys.argv[4]]


for file in vtkfiles:
   
   #Read the vtk file
   reader = vtk.vtkStructuredPointsReader()
   reader.ReadAllScalarsOn()
   reader.ReadAllVectorsOn()
   reader.SetFileName(file)
   reader.Update()

   # Set default y and z values
   y = 0.5
   z = 0.5

    # Create the line source to use for the probe lines.
   line = vtk.vtkLineSource()
   line.SetPoint1(x0,y,z)
   line.SetPoint2(x1,y,z)
   line.SetResolution(Nx)

   # Move the line into place and create the probe filter.  For
   # vtkProbeFilter, the probe line is the input, and the underlying data
   # set is the source.
   probe = vtk.vtkProbeFilter()
   probe.SetInputConnection(line.GetOutputPort())
   probe.SetSourceData(reader.GetOutput())
   probe.Update()
   data=probe.GetOutput()

   #Extract variables from point data
   ptdata   = data.GetPointData()
   arrayid  = ptdata.SetActiveVectors("velocity")
   velocity = ptdata.GetArray(arrayid)
   arrayid  = ptdata.SetActiveScalars("density")
   density = ptdata.GetArray(arrayid)
   arrayid  = ptdata.SetActiveScalars("pressure")
   pressure = ptdata.GetArray(arrayid)

   # Creating numpy array for data
   numpy_data = np.zeros((Nx,4))

   for i in range(Nx):
      x  = data.GetPoint(i)[0]
      v  = velocity.GetTuple3(i)[0]
      d  = density.GetTuple1(i)
      p  = pressure.GetTuple1(i)
      numpy_data[i] = np.array([x,d,v,p])

   
   fname = file[:-4] 
   fig,axs = plt.subplots(1,3,figsize=(20,5))
   axs[0].plot(numpy_data[:,0],numpy_data[:,1])
   axs[0].set_xlabel('x')
   axs[0].set_ylabel('density')
   axs[1].plot(numpy_data[:,0],numpy_data[:,2])
   axs[1].set_xlabel('x')
   axs[1].set_ylabel('velocity')
   axs[2].plot(numpy_data[:,0],numpy_data[:,3])
   axs[2].set_xlabel('x')
   axs[2].set_ylabel('pressure')
   plt.savefig(f"{fname}.pdf",bbox_inches='tight')
   plt.close()

   
   np.savetxt(f"{fname}.txt",numpy_data)
