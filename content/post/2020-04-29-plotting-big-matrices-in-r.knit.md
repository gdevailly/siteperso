---
title: Plotting heatmaps from big matrices in R
author: Guillaume Devailly
date: '2020-04-29'
slug: Plotting-big-matrices-in-r
categories:
    - R
tags:
    - R
    - data visualisaton
header:
    caption: ''
image: ''
math: true
draft: false
---



In genomics, and in many other -omics or big data fields, we often try to 
[plot big matrices](https://gdevailly.netlify.app/post/how-to-shuffle-ties-after-a-sort/).
By big matrices, I mean matrices that have more columns and/or rows than pixels in your device.
For example, if we have a 50 column-matrix with one line per gene (~20,000 lines),
then there are likely more lines in
this matrix than pixels in your screen - 1080 vertical pixels in an HD screen 
(unless you're reading this in a fancy future of hyper high definition).

The problem with plotting matrix that have more lines than pixels in the screen is precisely that:
each pixel will have to represent several values from the matrix.
R's default behaviour in that case might not be optimal for what we intend to show.

## An example with numeric data 

Let's start with the generation of some toy data, to mimic what one can obtain with 
[epigenomics data](https://gdevailly.netlify.app/post/how-to-shuffle-ties-after-a-sort/).
I'm trying to generate a column-centered signal, which becomes stronger along the rows, while keeping some randomness. I'll leave the code for reproducibility reasons, but 
you don't need to understand this part to understand what comes next:






















