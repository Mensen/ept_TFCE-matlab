.. -*- mode: rst -*-

ept_TFCE-matlab
===============

This matlab toolbox is designed for the statistical analysis of already
pre-processed M-EEG data. The primary underlying method used is the
'threshold-free cluster-enhancement' technique followed by permutation to
calculate the significance of the differences. The data may either be channel
by time; channel by frequency; or channel by frequency by time.

Get more information
^^^^^^^^^^^^^^^^^^^^

For more information on the statistical method we recommend looking at the paper entitled `Advanced EEG analysis using threshold-free cluster-enhancement and non-parametric statistics <http://www.sciencedirect.com/science/article/pii/S1053811912010300>`_.

**Abstract**

    Advances in EEG signal analysis and its combination with other
    investigative techniques make appropriate statistical analysis of large EEG
    datasets a crucial issue. With an increasing number of available channels
    and samples, as well as more exploratory experimental designs, it has
    become necessary to develop a statistical process with a high level of
    statistical integrity, signal sensitivity which nonetheless produces
    results which are interpretable to the common user. Threshold-free
    cluster-enhancement has recently been proposed as a useful analysis tool
    for fMRI datasets. This approach essentially takes into account both a data
    point's statistical intensity and neighbourhood to transform the original
    signal into a more intuitive understanding of ‘real’ differences between
    groups or conditions. Here we adapt this approach to optimally deal with
    EEG datasets and use permutation-based statistics to build an efficient
    statistical analysis. Furthermore we compare the results with several other
    non-parametric and parametric approaches currently available using
    realistic simulated EEG signals. The proposed method is shown to be
    generally more sensitive to the variety of signal types common to EEG
    datasets without the need for any arbitrary adjusting of parameters.
    Moreover, a unique p-value is produced for each channel-sample pair such
    that specific questions can still be asked of the dataset while providing
    general information regarding the large-scale experimental effects. 


Get the latest code
^^^^^^^^^^^^^^^^^^^

In order to download the latest stable code using git, type the following to the terminal::

  git clone git://https://github.com/Mensen/ept_TFCE-matlab.git


Install the toolbox
^^^^^^^^^^^^^^^^^^^

To use the toolbox simply add the directory where you downloaded the toolbox to the Matlab path using ``setpath``


Mex Files
^^^^^^^^^

Increase the speed of the more computationally demanding operations, this toolbox uses mex files which are compiled files written in the c-language. Although the toolbox includes the pre-compiled mex file for both unix and Windows systems, we cannot garauntee full compatibility with your specific system setup. For this reason we include the original .c files which can be used to compile the code again taking into account the specific system you are running under. To mex the files, navigate the dependencies folders in the toolbox and in the matlab command window type::

  mex ept_mex_TFCE2D.c

and::

  mex ept_mex_TFCE3D.c
  mex ept_mex_TFCE.c

You may need to initially configure your system to use a specific compiler which can be achieved by typing the following in the matlab command window::

  mex -setup

Troubleshooting
^^^^^^^^^^^^^^^

Use the "issues" section on GitHub to ask any questions relating to the usage of the toolbox that is not otherwise described in the documentation.
