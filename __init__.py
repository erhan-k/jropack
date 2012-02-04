"""Utility and application modules for jro data processing

jroread.py	---jro data file and header reading utilities
isr3Bb.py	---isr3Bb processing methods (imports jroread)
((mstisr2.py	---processing methods for mstisr2 data sets - future development))

Typical usage:
1) place jropack folder and data folders in some working directory
2) start python as "python -pylab" in the working directory
3) import jropack.isr3Bb as jro
4) jro.runTest() to produce and save spectrograms or jro.runTest('rti') to produce rti's etc...

jro. TAB in ipython for module methods...

"""

__all__=['jroread','isr3Bb','mstisr2','jrobeam','3Bb_beam','mstisr_beam']
