**1. kEpsH**

The k-epsilon based eddy viscosity model developed by [Fujihiro Hamba (2017)](https://pubs.aip.org/aip/pof/article-abstract/29/2/025103/834119/History-effect-on-the-Reynolds-stress-in-turbulent?redirectedFrom=fulltext). This model includes history effect in a pipe swirling flow and it's Reynolds stresses were derived in cylindrical coordinates. The model has been implemented in OpenFOAM software and successfully tested on Steenbergen’s experimental data.

**2. Motivation**

Swilirng flow is a complex flow encountered mainly in cylindrical flow domains. Therefore, it is naturally connected to physical analysis in cylindrical coordinates that results in a turbulence model in cylindrical coordinates. Since OpenFOAM uses Cartesian coordinates, it was challenging to find a solution of connecting turbulence model equations in cylindrical coordinates with OpenFoam solver in Cartesian coordinates, in order to ensure better predictions of pipe swirling flow.

**3. Target platform**

The code of this turbulence model is fully compatible with foam-extend-4.0 version of OpenFOAM software, and it has been rigorously tested ensuring good stability and convergence during numerical iterations.

**4. How to set model**
   
1 - Download the source code.

2 - Copy directories <code>kEpsH</code> and <code>Make</code> into directory on path

<code>~/foam/foam-extend-4.0/src/turbulenceModels/insompresible/RAS</code>

and replace existing files. On this path run command

<code>wmake libso</code>

to compile added model.

3 - Copy file <code>TensorTemplateI.H</code> into directory on path

<code>~/foam/foam-extend-4.0/src/primitives/Tensor</code>

and replace existing file.

4 - Download the mesh of test case from the link
   
   https://drive.google.com/drive/folders/1IPPsKB14MyTl_AMfqg7NkYmSOPBAHDR9?usp=drive_link
   
   and paste downloaded <code>polyMesh</code> directory into <code>constant</code> subdirectory of <code>testCase</code> directory.

5 - Run the test case.

**5. Results**

The model is tested on Steenbergen’s test pipe whose main part is presented in Fig. 1. The pipe is 21m long and has an internal diameter of D=2R=70mm. The swirl generator contains guide vanes that generate water swirling flow in the pipe. The central nozzle of the swirl generator prevents the appearance of reverse flow in the coaxial part of the pipe. Steenbergen has performed measurements in several cross sections downstream of the swirl generator for several flow regimes. The concentrated vortex swirling flow with inlet swirl strength of S<sub>0</sub>\=0.102 and Reynolds number Re<sub>D</sub>\=300000 is chosen for test of this model.

![Alt text](https://github.com/djordjenov/kEpsH/blob/main/testCase/TestRig.png)

Fig. 1. Test pipe with swirling genarator and positions of control cross sections.

This specific model intended for swirling flow simulations gives significantly better results compared to existing models derived in Cartesian coordinates in terms of axial and circumferential velocity prediction (see Fig. 2.). This model also gives better swirl decay predictions.

![Alt text](https://github.com/djordjenov/kEpsH/blob/main/testCase/VelocityProfiles.png)

Fig. 2. Axial and circumferential velocity profiles.
