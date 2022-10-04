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
