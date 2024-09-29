# Single Line Alignment with Context (SLAC)

SLAC is a simple Python tool built around an encoding system to enable visually intuitive, text-based representations of 3 aligned sequences, within a single continuous string of characters.

These 3 aligned sequences are intended to be:
- Query genomic sequence
- Query coding sequence
- Full hit relative to genomic sequence

In the context of sequence search workflows, this enables something like a query's BLAST hit on a non-reference genome to be shown with graphical context about the hit's location and polymorphisms relative to the query's coding and non-coding structure. 

A key feature of SLAC is the ability to generate abbreviated previews (miniSLAC) which allow it to generate 'thumbnail previews', suitable for embedding within common short-text-based environments such as spreadsheets and interactive plot hover text. This may greatly aid a user's ability to quickly navigate results within a unified interface.

This repo contains the core Python module for SLAC and is freely provided for use in the bioinformatics community. 

SLAC is available as an interactive website where aligned sequences can be provided directly. It also provides example sequences to better understand its use.
https://ryanzujic.github.io/slacweb/

### *Note: SLAC is currently in pre-release and a packaged release to PyPI will be available in the coming weeks. In the meantime, I'd encourage anyone interested to follow this repo and have a play with the interactive site to get a feel for SLAC.*

Anyone interested in using SLAC for their work or aiding its development is very welcome to reach out to chat.
