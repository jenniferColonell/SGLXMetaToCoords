# SGLXMetaToCoords
SpikeGLX metadata readers

Code to read SpikeGLX metadata for a given binary, and create a geometry map for spike sorting or other analysis.
The output format is set in the code -- see the code comments for details.

This code uses the ~snsGeomMap, present in data taken with SpikeGLX 20230202-phase30 or later. For older data, it reads the ~snsShankMap, and has hardcoded physical dimensions for all Neuropixels probe part numbers. There is an option to create a new metadata file appending the ~snsGeometryMap.


