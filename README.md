<img width="1904" height="1606" alt="pipeline" src="https://github.com/user-attachments/assets/d0c71717-3e7a-4393-b3ab-4324e88e04cc" />

### Overview
This deep mutational scanning pipeline couples fluorescence-activated cell sorting 
with next-generation sequencing to quantify the functional effects of protein variants at scale. 
Starting from a pooled library of yeast cells expressing RHO variants, the pipeline produces 
variant effect scores across multiple protein properties, including activity under different 
conditions and per-cell protein abundance.

### Quality Control
Raw sequencing reads are processed through a series of quality control steps to ensure only 
high-confidence variant calls are retained. This includes adapter trimming to remove library 
preparation artifacts, merging of Paired-End reads to reconstruct full-length amplicons, and 
filtering based on read quality scores and variant abundances to remove low-confidence or 
underrepresented sequences.

### Downsampling
To correct for biases introduced by unequal cell numbers across FACS bins, merged reads for 
each bin are downsampled to match the actual frequency of cells sorted into that bin. This 
normalization step ensures that variant scores reflect true biological differences in protein 
function rather than sorting artifacts. To account for stochastic variation introduced by 
downsampling itself, the process is repeated 30 times per condition and per protein property, 
generating a distribution of scores that captures sampling uncertainty.

### Variant Effect Scoring
For each downsample, variants are called and their separation across FACS bins is modelled 
to compute variant effect scores. This approach quantifies how strongly each variant shifts 
the fluorescence distribution relative to the wild type, providing a continuous measure of 
functional effect for each protein property tested. Scores from all 30 downsamples are then 
averaged to produce robust, combined variant effect scores, which are assembled into final 
variant effect maps.

### Analysis
The resulting variant effect maps enable comparative analysis of any two protein properties 
simultaneously — for example, plotting dark activation scores against light activation scores 
to identify variants that selectively alter one functional property while leaving another 
intact. Nonsense and missense variants can be distinguished in this space, with missense 
variants of interest clustering away from the loss-of-function nonsense baseline.
