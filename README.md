# PyVAMR
A Python package for visualizing animal mitochondrial rearrangements.

## Installation

```pip install git+https://github.com/thecgs/PyVAMR.git```

## Usage

### Example 1:

You can easily draw a circular mitogenome map in human, just like this, using the [OGDRAW](https://chlorobox.mpimp-golm.mpg.de/OGDraw.html) theme.

```python
import pyvamr

fig, ax = pyvamr.draw_circos_MT("NC_012920.1", colors="OGDRAW")
```

The result is shown in the figure below:

![Schema of quickprot algorithm](doc/Fig.1.png#pic_center)

<center>Fig.1 A circular mitogenome map in human</center>

### Parameter Details:

```python
help(pyvamr.draw_circos_MT)

Help on function draw_circos_MT in module pyvamr.drawMT:

draw_circos_MT(file, output=None, abbr=False, isfilename2species=False, colors='mitofish', radius=25, show_gene_label=True, gene_label_fontsize=7, gene_label_inner=False, show_info=True, info_fontsize=15, show_legend=True, legend_size=6, legend_postion=(1, -0.1), show_GC_circos=True, GC_circos_height=0.3, GC_circos_color='grey', GC_circos_bin=50, GC_circos_step=50, start='tRNA-Phe', axes=None, direction=-1, figsize=(10, 10), tidyname=False, add_id=False, dpi=300)
    Descripton:
        Draw a mitochondrial circos.
        
    Parameters：
        file: {str} a genbankfile or NCBI accession ID.
        output: {str} a path of figure save.
        abbr: {bool} whether to abbreviate species names.
        isfilename2species: {bool} whether filename convert to species.
        colors: {str, dict} themes such as, Chen, Tan, ogdraw, mitofish,
                            mitofish1, mitoz,  gggenes, chloroplot, grey, igv.
        radius: {int} radius of the circle.
        show_gene_label: {bool} show gene label.
        gene_label_fontsize: {int} gene label fontsize.
        gene_label_inner: {bool} whether the control gene label is inner.
        show_info: {bool} show species name, GC context, genome length, and gene count information.
        info_fontsize: {int} information fontsize.
        show_legend: {bool} show fig legend.
        legend_size: {int} legend size.
        legend_postion: {tuple} such as (x, y), x=0-1, y=0-1.
        show_GC_circos: {bool} show GC circos.
        GC_circos_height {0-1} height of GC_circos.
        GC_circos_color {color} color of GC_circos.
        GC_circos_step {int} step of GC_circos (bp).
        GC_circos_bin {int} bin of GC_circos (bp).
        start: {str} initial feature, such as, ND1, ND2, ND3, ND4, ND4L, ND5, ND6,
                     COX1, COX2, COX3, ATPase6, ATPase8, Cytb, tRNA-His, tRNA-Pro,
                     tRNA-Thr, tRNA-Trp, tRNA-Met, tRNA-Asp, tRNA-Ala, tRNA-Gln,
                     tRNA-Ile, tRNA-Arg, tRNA-Tyr, tRNA-Phe, tRNA-Lys, tRNA-Gly,
                     tRNA-Asn, tRNA-Leu, tRNA-Glu, tRNA-Val, tRNA-Cys, tRNA-Ser,
                     12S rRNA, 16S rRNA, D-loop.
        axes: {matplotlib.projections.polar.PolarAxes}
        figsize: {tuple} fig of size.
        direction: {1, -1} clockwise: -1, anticlockwise: 1
        tidyname: {bool} tidy gene name.
        add_id: {bool} Species add to accession id from NCBI.
        dpi: {int} dpi value. the resolution in dots per inch.
```



### Example 2:

According to a study by [Wang et al.](https://link.springer.com/article/10.1186/s40850-025-00239-x), compared with the ancestor of Stylommatophora, the mitochondrial genes of *M. pictum* exhibited multiple rearrangement events, while the mitochondrial genes of *S. arundinetorum* showed only minor differences. You can draw the two species separately for comparison (Fig. 2 A–B), or you can draw them together on the same figure, and using the [MitoFish](https://mitofish.aori.u-tokyo.ac.jp/annotation/draw) theme.

```python
import pyvamr
import  matplotlib.pyplot as plt

fig, axs = plt.subplots(1, 3, figsize=(20, 20/3), subplot_kw={'projection':'polar'})
plt.subplots_adjust(wspace=0.3)

pyvamr.draw_circos_MT(file="OP311642", colors="MitoFish", radius=20, gene_label_fontsize=5, gene_label_inner=False, show_info=True, show_legend=False, axes=axs[0], info_fontsize=6)
pyvamr.draw_circos_MT(file="OP289102", colors="MitoFish", radius=20, gene_label_fontsize=5, gene_label_inner=False, show_info=True, show_legend=False, axes=axs[1], info_fontsize=6)
pyvamr.draw_circos_MT(file="OP311642", colors="MitoFish", radius=20, gene_label_fontsize=5, gene_label_inner=False, show_info=False,show_legend=False, show_GC_circos=False, axes=axs[2])
pyvamr.draw_circos_MT(file="OP289102", colors="MitoFish", radius=12, gene_label_fontsize=4, show_gene_label=True,   show_info=False,show_legend=False, show_GC_circos=False, axes=axs[2])

pyvamr.add_tag(axs=axs, by_row=True)

axs[2].text(0.5, 0.5, s="Outer circle: Meghimatium pictum\nInner circle: Succinea arundinetorum", size=6, ha='center', va='center', style='italic')

fig.savefig("Fig.svg", bbox_inches='tight')
```

The result is shown in the figure below:

![Schema of quickprot algorithm](doc/Fig.2.png#pic_center)

<center>Fig.2 Three circular mitogenome map in two Stylommatophora species</center>

### Example 3:

Comparing the mitochondrial genomes of multiple species, Circos plots are clearly not the best choice; a study found that prevalent intraspecific gene rearrangements in [*Phrynocephalus*](https://www.mdpi.com/2073-4425/13/2/203). Below, we use a linear plot to visualize gene rearrangements. Use the gggenes theme and and open the circular mitochondrial genome at the specified site (such as, 12S rRNA).

```python
import pyvamr

accessions = ["KJ551842", "MF039062", "MF039063", "MF039065", "KJ630904", "MF039061",
              "MF039059", "MF039060", "MF039058", "KP126516", "MK284225", "KP232959",
              "KJ885621", "KJ830752", "KP279760", "KM093859", "MF039064", "KJ749841",
              "KM093858", "OL493803", "OL493804", "KC578685", "KC119493", "MK284224"]

pyvamr.draw_linear_MT_nonproportional(accessions, start='12S rRNA', species_offset=-0.0001, add_id=True, colors='GGGENES')
```

The result is shown in the figure below:

![Schema of quickprot algorithm](doc/Fig.3.png#pic_center)

<center>Fig.3 Linear plot in Phrynocephalus species</center>

### Example 4：

As shown in Fig.3, the genes are not drawn to scale. Below is a scaled drawing of the mitochondrial structure using the Tan theme.

```python
import pyvamr

accessions = ["MF039060", "MF039058", "KP126516", "MK284225"]
pyvamr.draw_linear_MT(accessions, start='12S rRNA', colors='Tan', add_id=True)
```

The result is shown in the figure below:

![Schema of quickprot algorithm](doc/Fig.4.png#pic_center)

<center>Fig.4 Linear plot of proportional in Phrynocephalus species</center>

## Themes

### Built-in themes:

It comes with 10 built-in themes, including Chen, Tan, OGDRAW, MitoFish, MitoFish1, MitoZ, Chloroplot, Grey, IGV, gggenes.

```python
import pyvamr
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 5, figsize=(15, 5), subplot_kw={'projection':'polar'})
axes = [a for ax in axs for a in ax]
themes = ["Chen", "Tan", "OGDRAW", "MitoFish", "MitoFish1", "MitoZ", "Chloroplot", "Grey", "IGV", "gggenes"]

for theme, ax in zip(themes, axes):
    pyvamr.draw_circos_MT("NC_012920.1", colors=theme, radius=12, show_gene_label=False, 
                              show_info=False, show_legend=False, show_GC_circos=False, axes=ax)
    ax.text(0.5, 0.5, s=theme, ha='center', va='center')
```

![Schema of quickprot algorithm](doc/Fig.5.png#pic_center)

<center>Fig.5 10 built-in themes</center>

### Custom colors:

You can also customize the color scheme to your liking. like this,

```python
MTColors_by_Set3 = {'source':"#000000", 'ND1': '#8DD3C7', 'ND2': '#8DD3C7', 'ND3': '#8DD3C7', 'ND4L':'#8DD3C7',
                    'ND4': '#8DD3C7', 'ND5': '#8DD3C7', 'ND6': '#8DD3C7', 'COX1':'#BEBADA', 'COX2':'#BEBADA',
                    'COX3':'#BEBADA', 'ATPase6':'#FDB462', 'ATPase8':'#FDB462', 'Cytb':'#80B1D3', 
                    'tRNA-His':'#FFFFB3', 'tRNA-Pro':'#FFFFB3', 'tRNA-Thr':'#FFFFB3', 'tRNA-Trp':'#FFFFB3',
                    'tRNA-Met':'#FFFFB3', 'tRNA-Asp':'#FFFFB3', 'tRNA-Ala':'#FFFFB3', 'tRNA-Gln':'#FFFFB3',
                    'tRNA-Ile':'#FFFFB3', 'tRNA-Arg':'#FFFFB3', 'tRNA-Tyr':'#FFFFB3', 'tRNA-Phe':'#FFFFB3',
                    'tRNA-Lys':'#FFFFB3', 'tRNA-Gly':'#FFFFB3', 'tRNA-Asn':'#FFFFB3', 'tRNA-Leu':'#FFFFB3',
                    'tRNA-Glu':'#FFFFB3', 'tRNA-Val':'#FFFFB3', 'tRNA-Cys':'#FFFFB3', 'tRNA-Ser':'#FFFFB3',
                    '12S rRNA':'#FB8072', '16S rRNA':'#FB8072', 'D-loop': '#E5C494'
                   }
```

