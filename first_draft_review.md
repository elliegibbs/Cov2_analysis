# Review of Analysis and Code repos

**date**: 2021-07-19

NOTE: No analyses as of 7/19

## Review Notes

- See instructions in Pull Request for incorporating changes
- Comments in julia code are preceded with `#`
- Comments in markdown look like this: `<!-- this is a comment -->`
- Comments that require action begin with `TODO`, others are just for information
- You may delete comments inside other files once you've addressed them
- Please keep this file in the repository - you may add your own responses if aplicable.

## Analysis plan

- Not in markdown
- Just includes functions, not a description of intentions, nor visuals

## Code repo

https://github.com/elliegibbs/FinalBisc195.jl

- No functions for analysis
  - New functions will need docstrings and tests
- Most functions from earlier assignments don't have docstrings

## Analysis repo

- Code doesn't run because of errors in Fasta parsing
- Your protein translation just starts from the beginning instead of finding ORF
  - I suggest you reconsider and download protein sequences instead
- Plots 