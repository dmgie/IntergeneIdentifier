# How to use it

To use this program, simply clone the entire repository , and place your reads within the same folder. Or alternatively, place the contents of your folder onto the parent folder of your reads i.e 

```
/
|script1.sh
|workflowt_script.sh
|----|reads
```



# Requirements

- Python3
- Bash
- Samtools
- R


# Notes:
It currently only has been tested on Linux (Ubuntu Server to be specific), so further testing is required. The `depth-add-name` script especially has been compiled for linux, so would require re-compiling for other. The Rust script used to compile the `depth-add-name` binary is in the repo: `dmgie/intergene-finder`
The seed_and_extend_IGR_folding.py script requires the exe files of the ViennaRNA package in a specific folder. If you want to use this script you might have to change the path to the exe file in the code. In the workflow.sh the normal seed_and_extend_IGR.py is called, which does not require the ViennaRNA package.
