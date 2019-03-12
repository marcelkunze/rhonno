# The Neural Network Objects

Johannes Steffens, Marcel Kunze, Helmut Schmücker

Neural Network Objects (NNO) is a C++ class library that implements the
most popular conventional neural networks together with novel
[incremental
models](http://www.ki.inf.tu-dresden.de/~fritzke/research/incremental.html)
that have been invented at Bochum University. The package is publicly
available and has proven versatile in a broad range of applications over
the past years. Recently NNO has been completely revised in order to
take full advantage of the [ROOT framework](http://root.cern.ch) for
data management and graphics. The Rho Analysis Framework is distributed with
all the example programs and training files mentioned in this document.


Architecture
============

At the time being the package comprises

| **Supervised Training Models**    | Multi-Layer Perceptron (TMLP,     |
|                                   | TXMLP)                            |
|                                   |                                   |
|                                   | Fisher Discriminant (TFD)\        |
|                                   | Supervised Growing Cell Structure |
|                                   | (TSGCS)\                          |
|                                   | Supervised Growing Neural Gas     |
|                                   | (TSGNG)                           |
|                                   |                                   |
|                                   | Neural Network Kernel (TNNK)      |
|-----------------------------------|-----------------------------------|
| **Unsupervised Training Models**  | Learning Vector Quantisation      |
|                                   | (TLVQ)\                           |
|                                   | Growing Cell Structure (TGCS)\    |
|                                   | Growing Neural Gas (TGNG)         |

The design foresees that all models are derived from the same abstract
base class ***VNeuralNet***. The common base class enforces a unique
interface to data management, training and recall cycles and graphics
operations at one central place. ***VSupervisedNet*** and
***VUnsupervisedNet*** both inherit from ***VNeuralNet*** and take care
of the different learning paradigms. In addition specific
implementations of the networks can utilize a plotter to produce a live
updating graphics window to control the training progress: The abstract
***VNeuralNetPlotter*** interface allows to plug in a graphics engine,
like the default **TSimpleNeuralNetPlotter**.

Implementation
==============

The ***VNeuralNet*** abstract interface defines the following contract
for the implementation of specific neural network models:

  ------------------------------------------------------------------------
  ## Abstract interface for all networks
  
  virtual void AllocNet() = 0;
  
  virtual void InitNet() = 0;
  
  virtual void WriteText() = 0;
  
  virtual void WriteBinary() = 0;
  
  virtual void ReadText() = 0;
  
  virtual void ReadBinary() = 0;
  
  virtual Double\_t\* Recall(NNO\_INTYPE\* in,NNO\_OUTTYPE\* out=0) = 0;
  
  virtual Double\_t Train(NNO\_INTYPE\* in,NNO\_OUTTYPE\* out=0) = 0;
  
  ------------------------------------------------------------------------

*AllocNet* acquires resources and is executed during network
construction. *InitNet* sets up the network weight matrix. *WriteText*
persists the network as an ASCII file, *WriteBinary* produces a binary
file. The corresponding reading versions are able to regenerate a
network in the state it was at the time when it was saved. The *Recall*
method takes an input vector as parameter and returns the corresponding
network output. The *Train* function takes a pair of input/output
vectors, performs a *Recall*, modifies the weight matrix to better adapt
the input probability density function and returns the squared error of
the sample.

Besides the abstract interface, concrete methods have been implemented
to support the execution of training cycles and to set training
parameters:

  -------------------------------------------------------------
  
  ## Training and testing
  
  Double\_t TrainEpoch(TDataServe \*server, Int\_t nEpoch=1);
  
  Double\_t TestEpoch(TDataServe \*server);
  
  void BalanceSamples(Bool\_t yesNo = kTRUE);
  
  virtual void SetMomentumTerm(Double\_t f);
  
  virtual void SetFlatSpotElimination(Double\_t f);
  
  -------------------------------------------------------------

**TDataServe** is a mini database to support management of input/output
vector relations. It allows to partition datasets into training and test
samples, retrieve arbitrary samples and shuffle a data set prior to a
new training cycle. *TrainEpoch* and *TestEpoch* are functions to train
and test networks with a complete set of vectors out of a **TDataServe**
object. The *BalanceSamples* option allows to have equal training
statistics for good and bad samples, independent of the number of
vectors per sample. The application of a momentum term might lead to
faster convergence in some applications by noting the direction of
gradient descent, the flat spot elimination might improve training
progress in regions where the derivatives of the error matrix are near
zero.

Each network implementation has to implement the abstract interface
mentioned above. As an example for the integration of an independent
neural network implementation into the context of NNO we have managed to
support J.P. Ernenwein’s Neural Network Kernel: The ***TNNK*** interface
yields seamless access to the Neural Network Kernel in the scope of NNO.

NetworkTrainer
==============

Network training requires to identify pairs of input vectors and output
vectors out of a dataset of good and bad samples to describe the problem
at hand. It usually takes a certain amount of time to select suiting
quantities and assemble corresponding training and test files prior to
network training and write a corresponding training program or macro.
However, it turns out that ROOT is performing enough to allow for
interactive training of large networks with large data samples out of
arbitrary ROOT files in one go. In that spirit a *NetworkTrainer*
program has been written on the basis of the NNO package.
*NetworkTrainer* assists to

-   Assemble training and testing data sets out of ROOT trees

-   Define the network architecture

-   Define a training schedule

-   Persist networks

-   Generate C++ code to perform network recall

At the time being *NetworkTrainer* reads an ASCII steering file when it
launches (a GUI is in preparation). The steering file knows the
following directives:

---------------------------------------------------------

+-----------------------+-----------------------+-----------------------+
| ### *Parameter*       | ***Type***            | #### Description      |
|                       |                       |                       |
|                       | ***I = input\         |                       |
|                       | O = output***         |                       |
|                       |                       |                       |
|                       | ***H = hidden***      |                       |
|                       |                       |                       |
|                       | ***C = cells***       |                       |
+-----------------------+-----------------------+-----------------------+
| ***fisher***          | vector\               | Multi-layer           |
|                       | (I O)                 | perceptron (0 hidden  |
|                       |                       | layer)                |
+-----------------------+-----------------------+-----------------------+
| ***mlp***             | vector\               | Multi-layer           |
|                       | (I H O)               | perceptron (1 hidden  |
|                       |                       | layer)                |
+-----------------------+-----------------------+-----------------------+
| ***xmlp***            | vector\               | Multi-layer           |
|                       | (I H H O)             | perceptron (2 hidden  |
|                       |                       | layers)               |
+-----------------------+-----------------------+-----------------------+
| ***tnnk***            | vector\               | Multi-layer           |
|                       | (I H H O)             | perceptron (Neural    |
|                       |                       | Network Kernel)       |
+-----------------------+-----------------------+-----------------------+
| ***sgng***            | vector\               | Supervised growing    |
|                       | (I C O)               | neural gas            |
+-----------------------+-----------------------+-----------------------+
| ***sgcs***            | vector\               | Supervised growing    |
|                       | (I C O)               | cell structures       |
+-----------------------+-----------------------+-----------------------+
| ***gng***             | vector\               | Growing neural gas    |
|                       | (I C)                 |                       |
+-----------------------+-----------------------+-----------------------+
| ***gcs***             | vector\               | Growing cell          |
|                       | (I C)                 | structures            |
+-----------------------+-----------------------+-----------------------+
| ***lvq***             | vector\               | Learning vector       |
|                       | (I C)                 | quantization          |
+-----------------------+-----------------------+-----------------------+
| ***start***           | int                   | First training epoch  |
+-----------------------+-----------------------+-----------------------+
| ***stop***            | int                   | Last training epoch   |
+-----------------------+-----------------------+-----------------------+
| ***epoch***           | int                   | Number of training    |
|                       |                       | samples per epoch     |
+-----------------------+-----------------------+-----------------------+
| ***test***            | int                   | Number of test        |
|                       |                       | samples per epoch     |
+-----------------------+-----------------------+-----------------------+
| ***networkpath***     | string                | Directory to save the |
|                       |                       | trained networks      |
+-----------------------+-----------------------+-----------------------+
| ***datapath***        | string                | Directory to look up  |
|                       |                       | data files            |
+-----------------------+-----------------------+-----------------------+
| ***file***            | string                | ROOT training file    |
|                       |                       | containing good and   |
|                       |                       | bad samples           |
+-----------------------+-----------------------+-----------------------+
| ***pro***             | string                | ROOT training file    |
|                       |                       | containing good       |
|                       |                       | samples (1D output    |
|                       |                       | only)                 |
+-----------------------+-----------------------+-----------------------+
| ***con***             | string                | ROOT training file    |
|                       |                       | containing bad        |
|                       |                       | samples (1D output    |
|                       |                       | only)                 |
+-----------------------+-----------------------+-----------------------+
| ***tree***            | string                | ROOT tree that acts   |
|                       |                       | as source to assemble |
|                       |                       | the vectors           |
+-----------------------+-----------------------+-----------------------+
| ***cut***             | string                | ROOT TFormula for     |
|                       |                       | preselection of       |
|                       |                       | samples               |
+-----------------------+-----------------------+-----------------------+
| ***input***           | string                | Input vector, ROOT    |
|                       |                       | TFormulae (separated  |
|                       |                       | by colon)             |
+-----------------------+-----------------------+-----------------------+
| ***output***          | string                | Output vector, ROOT   |
|                       |                       | TFormulae (separated  |
|                       |                       | by colon)             |
+-----------------------+-----------------------+-----------------------+
| ***transfer***        | string                | Transfer function     |
|                       |                       | (TR\_FERMI,TR\_LINEAR |
|                       |                       | ,TR\_LINEAR\_BEND,TR\ |
|                       |                       | _SIGMOID)             |
+-----------------------+-----------------------+-----------------------+
| ***momentum***        | float                 | Momentum term         |
+-----------------------+-----------------------+-----------------------+
| ***scale***           | float                 | Global scale factor   |
|                       |                       | to apply to input     |
|                       |                       | layer                 |
+-----------------------+-----------------------+-----------------------+
| ***inscale***         | vector                | Scale factors to      |
|                       |                       | apply to input layer  |
+-----------------------+-----------------------+-----------------------+
| ***outscale***        | vector                | Scale factors to      |
|                       |                       | apply to output layer |
+-----------------------+-----------------------+-----------------------+
| ***autoscale***       | bool                  | Determine scale       |
|                       |                       | factors to apply to   |
|                       |                       | input layer           |
+-----------------------+-----------------------+-----------------------+
| ***plot***            | bool                  | Produce graphics      |
|                       |                       | output (1D output     |
|                       |                       | only)                 |
+-----------------------+-----------------------+-----------------------+
| ***balance***         | bool                  | Enforce presentation  |
|                       |                       | of equal number of    |
|                       |                       | good and bad samples  |
+-----------------------+-----------------------+-----------------------+

---------------------------------------------------------

A sample steering file for training of a selector to separate different
charged particles in a typical HEP experiment could look like the
following:

  ---------------------------------------------------------
  
  \# Training of PIDSelectors with NNO
  
  \#define the network topology
  
  xmlp 7 15 10 1
  
  transfer TR\_FERMI
  
  momentum 0.2
  
  balance true
  
  plots true
  
  test 10000
  
  start 1
  
  stop 200
  
  \#define the data source
  
  datapath ../Data
  
  networkpath ../Networks
  
  file PidTuple1.root
  
  file PidTuple2.root
  
  \#set up the input layer (use branch names)
  
  tree PidTuple
  
  cut mom&gt;0.5&&dch&gt;0&&dch&lt;10000
  
  input mom:acos(theta):svt:emc:drc:dch:ifr:ifrExp:ifrAdd
  
  autoscale true
  
  \#set up the output layer (use branch names)
  
  \#Particles pid = {electron=1,muon,pion,kaon,proton}
  
  output abs(pid)==3
  
  \#end of file
  
  ---------------------------------------------------------

The example above reads two input files, assembles a data server using
all samples surviving the cut and runs for 200 training epochs with a
7-15-10-1 multi-layer perceptron using all available samples. In the
course of the training after each epoch a persistent network file
NNOxxxx.TXMLP is saved into the Networks directory, where xxxx denotes
the epoch number. At the end, NetworkTrainer produces a template recall
function that can be plugged into another program that wants to make use
of a network. For the above example the file RecallTXMLP.cpp looks like
is shown below for illustration purposes:

  -----------------------------------------------------------------------
  
  // TXMLP network trained with NNO NetworkTrainer at Fri Apr 27
  
  // Input parameters mom:acos(theta):svt:emc:drc:dch:ifr:ifrExp:ifrAdd
  
  // Output parameters abs(pid)==3
  
  // Training files:
  
  //../Data/PidTuple1.root
  
  //../Data/PidTuple2.root
  
  \#include "RhoNNO/TXMLP.h"
  
  Double\_t\* Recall(Double\_t \*invec)
  
  {
  
  static TXMLP net("TXMLP.net");
  
  Float\_t x\[7\];
  
  x\[0\] = 0.76594 \* invec\[0\]; // mom
  
  x\[1\] = 2.21056 \* invec\[1\]; // acos(theta)
  
  x\[2\] = 0.20365 \* invec\[2\]; // svt
  
  x\[3\] = 2.2859 \* invec\[3\]; // emc
  
  x\[4\] = 1.75435 \* invec\[4\]; // drc
  
  x\[5\] = 0.00165 \* invec\[5\]; // dch
  
  x\[6\] = 0.85728 \* invec\[6\]; // ifr
  
  return net.Recall(x);
  
  }
  
  -----------------------------------------------------------------------

How to get and build RhoNNO
===========================

Rho resides on github. If you wish to install and
build, you have to check out and build the corresponding packages and applications:

&gt; git clone https://github.com/marcelkunze/rhonno

&gt; cd rhonno/RhoNNO

&gt; make

&gt; ../bin/strain

&gt; ../bin/NetworkTrainer pid.nno

Prior to start you should install the most recent production version of
ROOT from root.cern.ch and set the ROOTSYS environment variable
correspondingly. In addition you have to add \$ROOTSYS/lib and \$RHO/lib
to your LD\_LIBRARY\_PATH in order to resolve the shared libs.

References
==========

[The Neural Network
Objects](https://github.com/marcelkunze/rhonno/blob/master/doc/RhoNNO.pdf), J.Steffens,
M.Kunze
