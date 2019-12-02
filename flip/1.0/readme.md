## FLIP (Fluorescence Imaging Pipeline)

### Getting Started 
---

#### CLI (command line interface)

To run the cli, run 

```
$ python FLIP.py -d <ps2 collection path> -o <output folder>
```

To run the docker container, type

```
$ docker run --rm -v $PWD:/data -w /data cyverse/flip/flip:1.0 -d <ps2 collection path> -o final_out
```

This will run through the binary to png conversion, run multithreshold image segmentation, and then generate aggregate and fluorescence files all in one go and move the output files into the output folder.
