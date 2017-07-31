# adjusted LHAPDF

Should the code not compile when including `LHAPDF`, need to use the patched version here. To use it do the following:

Edit file `$CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/lhapdf.xml` to adjust the path of the *include* files of the LHAPDF installation in which this README resides (i.e. only change the path for `environment name="INCLUDE"`). You need to use the full path, and the file should look similar to the example below:
```
<tool name="lhapdf" version="6.1.6">
  <lib name="LHAPDF"/>
  <client>
    <environment name="LHAPDF_BASE" default="/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/lhapdf/6.1.6"/>
    <environment name="LIBDIR" default="$LHAPDF_BASE/lib"/>
    <environment name="INCLUDE" default="/afs/cern.ch/work/c/clange/MSSMbbH/CMSSW_8_0_20_patch1/src/Analysis/lhapdf/6.1.6/include"/>
  </client>
  <runtime name="LHAPDF_DATA_PATH" value="$LHAPDF_BASE/share/LHAPDF"/>
  <use name="yaml-cpp"/>
  <runtime name="ROOT_INCLUDE_PATH" value="$INCLUDE" type="path"/>
  <use name="root_cxxdefaults"/>
</tool>
```
Then run
```
scram setup lhapdf
```
Make sure there are no complaints/errors.

To check if things worked:
```
scram tool info lhapdf
```

Then recompile your code from the `Analysis` directory:
```
scram b clean && scram b -j8
```
