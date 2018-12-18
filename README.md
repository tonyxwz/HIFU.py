# pyHIFU

HIFU simulation in python, with ray tracing algorithm: Trident

## Goals

- prove for cubes and cylinders

## routine

```algorithm
Initialize the voxel coordinates of the whole system
Initialize the geometric representation of the object

for each element in trasducer:
    for each trident from the element:
        1. find closest intersection with geometry surfaces
            (loop through all the volumes, compare length and attach volume index)
        2. calculate all intersection with voxel faces -> list
            (from 1 we have the end of current ray, we should also have the start of the ray
             together, all the voxels on the path should be determined)
        3. accumulate the influence of current trident to all the voxels affected
```

## prerequisite modules

- geometric module to describe volumes

- physics module to calculate outgoing angles, reflection angles

- pyHIFU.HIFU: run the whole HIFU system

- pyHIFU.ray: US ray and trident ray

- pyHIFU.transducer: transducer matrix, elements

## next steps

finite element method mesh-like surfaces in geometric modules