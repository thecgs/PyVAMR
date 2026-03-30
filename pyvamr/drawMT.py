#!/usr/bin/env python
# coding: utf-8

import math, re
import logging
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyArrow
from matplotlib.lines import Line2D
from .parserGB import get_features
from .config import MTColors_legends, FullName2AbbrName

def rotation_text(theta):
    d =  (theta * 180 / math.pi)
    if d <180:
        return 90 - d 
    else:
        return  90 -d + 180

def get_GC(seq):
    return round((seq.count("G") + seq.count("C")) / (seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C"))*100, 2)     

def get_GC_bar_param(features, bin=50, step=50):
    x = []
    y = []
    w = []
    scale_factor = (2 * math.pi) / (len(features[0].mtgenome) + 1)
    for i in range(0, len(features[0].mtgenome), step):
        if (i + bin) < len(features[0].mtgenome):
            start = i + 1
            end = i + bin + 1
            theta = ((end - start) / 2 + start) * scale_factor
            width = (end - start) * scale_factor
            GC = get_GC(features[0].mtgenome[i: i+bin])
            y.append(GC/100)
            x.append(theta)
            w.append(width)
        else:
            start = i + 1
            end = len(features[0].mtgenome) + 1
            theta = ((end - start) / 2 + start) * scale_factor
            width = (end - start) * scale_factor
            GC = get_GC(features[0].mtgenome[i: len(features[0].mtgenome)])
            y.append(GC/100)
            x.append(theta)
            w.append(width)
    return x, y, w
    
def stat_features(features):
    tRNA_count = 0
    rRNA_count = 0
    PCG_count = 0
    for feature in features:
        if re.search('tRNA', feature.name) and feature.join==False:
            tRNA_count += 1
        elif re.search('rRNA', feature.name) and feature.join==False:
            rRNA_count += 1
        elif re.search('COX|ATP|Cytb|ND', feature.name) and feature.join==False:
            PCG_count += 1

    for i in {i.name for i in features if i.join}:
        if re.search('tRNA', feature.name) and feature.join==False:
            tRNA_count += 1
        elif re.search('rRNA', feature.name) and feature.join==False:
            rRNA_count += 1
        elif re.search('COX|ATP|Cytb|ND', feature.name) and feature.join==False:
            PCG_count += 1
    return tRNA_count, rRNA_count, PCG_count

def draw_circos_MT(file,
                   output=None,
                   abbr=False,
                   isfilename2species=False,
                   colors="mitofish",
                   radius=25,
                   show_gene_label = True,
                   gene_label_fontsize=7,
                   gene_label_inner=False,
                   show_info=True,
                   info_fontsize=15,
                   show_legend=True,
                   legend_size=6,
                   legend_postion=(1, -0.1),
                   show_GC_circos=True,
                   GC_circos_height = 0.3,
                   GC_circos_color = 'grey',
                   GC_circos_bin = 50,
                   GC_circos_step = 50,
                   start="tRNA-Phe",
                   axes=None,
                   direction=-1,
                   figsize=(10,10),
                   tidyname=False,
                   add_id = False
                  ):
    """
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
    """
    
    if axes==None:
        fig, ax = plt.subplots(1,1, subplot_kw={'projection':'polar'}, figsize=figsize)
    else:
        ax = axes
    
    ax.set_theta_zero_location('N')
    #ax.set_theta_offset(offset)
    ax.set_theta_direction(direction)
    ax.grid(False)
    ax.set_frame_on(False)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
    features = get_features(file, abbr=abbr, isfilename2species=isfilename2species, colors=colors, start=start)
    
    try:
        #str(features[0].mtgenome)  # UndefinedSequenceError
        GC = get_GC(features[0].mtgenome)
        genome_length = len(features[0].mtgenome)
        tRNA_count, rRNA_count, PCG_count = stat_features(features)
        if add_id:
            info = features[0].accession +"; " + features[0].name+f"\n{genome_length} bp\nGC: {GC}%\n{PCG_count} PCGs; {rRNA_count} rRNAs; {tRNA_count} tRNAs"
        else:
            info = features[0].name+f"\n{genome_length} bp\nGC: {GC}%\n{PCG_count} PCGs; {rRNA_count} rRNAs; {tRNA_count} tRNAs"
    except:
        #show_info = False
        show_GC_circos = False
        logger = logging.getLogger(__name__) 
        logger.setLevel(logging.DEBUG)
        logger.warning(features[0].file+" Sequence content is undefined.")
        tRNA_count, rRNA_count, PCG_count = stat_features(features)
        if add_id:
            info = features[0].accession +"; " + features[0].name+f"\n{PCG_count} PCGs; {rRNA_count} rRNAs; {tRNA_count} tRNAs"
        else:
            info = features[0].name+f"\n{PCG_count} PCGs; {rRNA_count} rRNAs; {tRNA_count} tRNAs"
        
    if show_info:            
        #GC = get_GC(features[0].mtgenome)
        #genome_length = len(features[0].mtgenome)
        ax.text(0.5, 0.5, s=info,
                fontsize=info_fontsize, style='italic', va='center', ha='center')
                
    if tidyname:
        for i, f in enumerate(features):
            features[i].name = FullName2AbbrName.get(f.name, f.name)
    
    for feature in features[1:]:
        scale_factor = (2 * math.pi) / (int(features[0].location.end) + 1)
        start = int(feature.location.start) + 1
        end = int(feature.location.end) + 1
        theta = ((end - start) / 2 + start) * scale_factor
        width = (end - start) * scale_factor
        
        if gene_label_inner==False:
            if feature.location.strand == 1:
                ax.bar(theta, 1, width=width, bottom=radius-5, color=feature.color, label=feature.name, linewidth=0.5, edgecolor='black')
                if show_gene_label:
                    ax.annotate(feature.name,  xy=(theta, radius-4), xytext=(theta, radius-1), arrowprops=dict(arrowstyle="-", connectionstyle="arc3"), 
                                horizontalalignment='center', verticalalignment='center', fontsize=gene_label_fontsize, rotation=rotation_text(theta))
            else:
                ax.bar(theta, 1, width=width, bottom=radius-6, color=feature.color, label=feature.name, linewidth=0.5, edgecolor='black')
                if show_gene_label:
                    ax.annotate(feature.name,  xy=(theta, radius-5), xytext=(theta, radius), arrowprops=dict(arrowstyle="-", connectionstyle="arc3"), 
                                horizontalalignment='center', verticalalignment='center', fontsize=gene_label_fontsize, rotation=rotation_text(theta))
            ax.plot([i* math.pi/180 for i in range(0, 361)], [radius-5]*361, c='black')
        
        else:
            if feature.location.strand == 1:
                ax.bar(theta, 1, width=width, bottom=radius, color=feature.color, label=feature.name, linewidth=0.5, edgecolor='black')
                if show_gene_label:
                    ax.annotate(feature.name,  xy=(theta, radius), xytext=(theta, radius-4), arrowprops=dict(arrowstyle="-", connectionstyle="arc3"), 
                                horizontalalignment='center', verticalalignment='center', fontsize=gene_label_fontsize, rotation=rotation_text(theta))
            else:
                ax.bar(theta, 1, width=width, bottom=radius-1, color=feature.color, label=feature.name, linewidth=0.5, edgecolor='black')
                if show_gene_label:
                    ax.annotate(feature.name,  xy=(theta, radius-1), xytext=(theta, radius-5), arrowprops=dict(arrowstyle="-", connectionstyle="arc3"), 
                                horizontalalignment='center', verticalalignment='center', fontsize=gene_label_fontsize, rotation=rotation_text(theta))
                
            ax.plot([i* math.pi/180 for i in range(0, 361)], [radius]*361, c='black')
            
    if colors == None:
        legend_elements = MTColors_legends['MITOFISH']
        ncol = 1
    elif isinstance(colors, str):
        legend_elements = MTColors_legends.get(colors.upper(), 'MITOFISH')
        if colors.upper() == "IGV":
            ncol = 4
        else:
            ncol = 1
    elif isinstance(colors, dict):
        legend_elements = [Patch(facecolor=colors[i], edgecolor='black', label=i) for i in colors if i!="source"]
        ncol = 4
    
    if show_legend:
        ax.legend(handles=legend_elements,
                  loc='center', 
                  bbox_to_anchor=legend_postion, 
                  ncol=ncol,
                  #frameon=False,
                  #shadow=False, 
                  prop={'size': legend_size},
                  title='')
    
    if show_GC_circos:
        x,y,w = get_GC_bar_param(features, step=GC_circos_step, bin=GC_circos_bin)
        ax.bar(x, [i*10 *GC_circos_height for i in y], width=w, bottom=radius-10, color=GC_circos_color, linewidth=0, edgecolor='black')
    
    if output!=None:
        fig.savefig(output, bbox_inches='tight')
    
    if axes==None:
        return fig, ax
 
def _get_mt_rect_ax(ax, features, cex=18000, show_info=True, info_fontsize=10, show_gene_label=True, gene_label_fontsize=5, add_id=False):
    #cex = 18000      # genome max length
    height = 0.2     # gene width
    linewidth = 0.03 # genome width
    
    for feature in features:
        start = int(feature.location.start)/cex
        end = int(feature.location.end)/cex
        
        if feature.type == 'source':
            rect = plt.Rectangle(xy=(0, 0.45), width=end, height=linewidth,  color=feature.color)
            species = feature.name
            genome_length = int(feature.location.end)
            
        elif feature.location.strand == 1:
            rect = plt.Rectangle(xy=(start, 0.45+linewidth), 
                                 width=end-start, 
                                 height=height, 
                                 facecolor=feature.color,
                                 edgecolor='black')
        else:
            rect = plt.Rectangle(xy=(start, 0.25), 
                                 width=end-start,
                                 height=height, 
                                 facecolor=feature.color,
                                 edgecolor='black')
                 
        ax.add_patch(rect)
        
        if feature.type != 'source' and show_gene_label:
            if feature.location.strand == 1:
                x = (end-start)/2 + start
                ax.text(x, y=0.45+height+linewidth+0.05, s=feature.name, rotation = 'vertical', size=gene_label_fontsize, ha='center')
            else:
                x = (end-start)/2 + start
                ax.text(x, y=0.45-height-linewidth+0.005, s=feature.name, rotation = 'vertical', size=gene_label_fontsize, va='top', ha='center')
                
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    if add_id:
        ax.text(-0.01, 0.5, s=species+" "+features[0].accession, ha='right', style='italic', fontsize=info_fontsize)
    else:
        ax.text(-0.01, 0.5, s=species, ha='right', style='italic', fontsize=info_fontsize)
    
    if show_info:
        ax.text(genome_length/cex+0.01, 0.55, s="+", ha='center', weight='bold', fontsize=info_fontsize)
        ax.text(genome_length/cex+0.01, 0.25, s="-", ha='center', weight='bold', fontsize=info_fontsize)
        ax.text(genome_length/cex+0.02, 0.47, s=str(genome_length)+'bp', ha='left', fontsize=info_fontsize)
    return None

def draw_linear_MT(files,
                   output=None,
                   abbr=False,
                   isfilename2species=False,
                   colors="mitofish",
                   show_gene_label = True,
                   gene_label_fontsize=5,
                   show_info=True,
                   info_fontsize=10,
                   show_legend=True,
                   legend_size=10,
                   legend_postion=(1, 0.5),
                   start="tRNA-Phe",
                   show_xaxis=True,
                   xaxisfontsize=10, 
                   hspace=0.5,
                   subplot_height_cex=1.5,
                   tidyname=False,
                   add_id = False
                  ):
    
    """
    Descripton:
        The order of mitochondrial genes, arranged in proportion to their size.
    
    Parameters：
        files: {str, list, tuple} one or more genbankfile or NCBI accession ID.
        output: {str} a path of fig save.
        abbr: {bool} whether to abbreviate species names.
        isfilename2species: {bool} whether filename convert to species.
        colors: {str, dict} themes such as, Chen, Tan, ogdraw, mitofish,
                            mitofish1, mitoz,  gggenes, chloroplot, grey, igv.

        show_gene_label: {bool} show gene label.
        gene_label_fontsize: {int} gene label fontsize.
        show_info: {bool} show genome length, and strand.
        info_fontsize: {int} information fontsize.

        show_legend: {bool} show fig legend.
        legend_size: {int} legend size.
        legend_postion: {tuple} such as (x, y), x=0-1, y=0-1.
        start: {str} initial feature, such as, ND1, ND2, ND3, ND4, ND4L, ND5, ND6,
                     COX1, COX2, COX3, ATPase6, ATPase8, Cytb, tRNA-His, tRNA-Pro,
                     tRNA-Thr, tRNA-Trp, tRNA-Met, tRNA-Asp, tRNA-Ala, tRNA-Gln,
                     tRNA-Ile, tRNA-Arg, tRNA-Tyr, tRNA-Phe, tRNA-Lys, tRNA-Gly,
                     tRNA-Asn, tRNA-Leu, tRNA-Glu, tRNA-Val, tRNA-Cys, tRNA-Ser,
                     12S rRNA, 16S rRNA, D-loop.
        show_xaxis: {bool} show x-axis.
        xaxisfontsize: {int} x-axis trick label fontsize.
        hspace: {0-1} spacing between subplots.
        subplot_height_cex: {float} a factors that control subplot scaling.
        tidyname: {bool} tidy gene name.
        add_id: {bool} Species add to accession id from NCBI.
    """    

    if not (isinstance(files, list)) and (not isinstance(files, tuple)):
        files = [files]
    
    genome_max_length = 0
    genomes = []
    for file in files:
        features = get_features(file, abbr=abbr, isfilename2species=isfilename2species, colors=colors, start=start)
        if tidyname:
            for i, f in enumerate(features):
                features[i].name = FullName2AbbrName.get(f.name, f.name)
        genomes.append(features)
        if genome_max_length < len(features[0].location):
            genome_max_length = len(features[0].location)
    
    if show_xaxis:
        fig, ax= plt.subplots(len(genomes)+1, 1, figsize=(20, len(genomes)*subplot_height_cex))
    else:
        fig, ax= plt.subplots(len(genomes), 1, figsize=(20, len(genomes)*subplot_height_cex))
        
    plt.subplots_adjust(hspace=hspace)
    
    if show_xaxis == False and len(genomes) == 1:
        for i, features in enumerate(genomes):
            _get_mt_rect_ax(ax=ax, features=features, cex=genome_max_length+1000, show_info=show_info,
                           info_fontsize=info_fontsize, show_gene_label=show_gene_label, 
                           gene_label_fontsize=gene_label_fontsize, add_id=add_id)
    else:
        for i, features in enumerate(genomes):
            _get_mt_rect_ax(ax=ax[i], features=features, cex=genome_max_length+1000, show_info=show_info,
                            info_fontsize=info_fontsize, show_gene_label=show_gene_label,
                            gene_label_fontsize=gene_label_fontsize, add_id=add_id)
        
    if show_xaxis:
        ax[-1].get_yaxis().set_visible(False)
        ax[-1].spines['left'].set_visible(False)
        ax[-1].spines['right'].set_visible(False)
        ax[-1].spines['bottom'].set_visible(False)
        ax[-1].set_xticks([i for i in range(0 , genome_max_length+1000, 1000)], labels=[i for i in range(0 , genome_max_length+1000, 1000)], **{"fontsize":xaxisfontsize})
        ax[-1].tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, direction='in', pad=-20)
    
    if colors == None:
        legend_elements = MTColors_legends['MITOFISH']
    elif isinstance(colors, str):
        legend_elements = MTColors_legends.get(colors.upper(), 'MITOFISH')

    if show_legend:
        ax[-1].legend(handles=legend_elements,
                      loc='upper right', 
                      bbox_to_anchor=legend_postion, 
                      ncol=12,
                      #frameon=False,
                      #shadow=False,
                      prop={'size': legend_size},
                      title='')
    
    if output!=None:
        fig.savefig(output, bbox_inches='tight')    
    return fig, ax
        
def get_box_param(genename):
    head_length=0.0025
    if genename in ["F","V","L","I","Q","M","W","A","N","C","Y","S","D","K","G","R","H","S","L","E","T","P"]:
        tail_length = 0.005
        text_offset = tail_length/2
        
    elif genename in ["ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"]:
        tail_length = 0.008
        if genename=="ND4L":
            text_offset = tail_length/8
        else:
            text_offset = tail_length/6
        
    elif genename in ["12S", "16S"]:
        tail_length = 0.005
        text_offset = tail_length/6
    elif genename in ['COX1','COX2','COX3']:
        tail_length = 0.008
        text_offset = tail_length/8
    elif genename in ["ATP6", "ATP8"]:
        tail_length = 0.008
        text_offset = tail_length/6
    elif genename == "NCR":
        tail_length = 0.007
        text_offset = tail_length/6
    elif genename == "Cytb":
        tail_length = 0.008
        text_offset = tail_length/5
    else:
        tail_length = 0.015
        text_offset = tail_length/4
    box_width = head_length + tail_length
    return head_length, tail_length, box_width, text_offset
    
def draw_linear_MT_nonproportional(files, output=None, abbr=False, isfilename2species=False, colors="mitofish", gene_label_fontsize=9, 
                                   gene_label_color="black", show_legend=True, legend_size=8, legend_postion=(1, 0),
                                   start="tRNA-Phe", height=0.6, species_fontsize=10, species_offset=-0.01, axes=None, add_id=False,
                                  ):
    """
    Descripton:
        The order of mitochondrial genes is not depicted in proportion to their gene size.
    
    Parameters：
        files: {str, list, tuple} one or more genbankfile or NCBI accession ID.
        output: {str} a path of fig save.
        abbr: {bool} whether to abbreviate species names.
        isfilename2species: {bool} whether filename convert to species.
        colors: {str, dict} themes such as, Chen, Tan, ogdraw, mitofish,
                            mitofish1, mitoz,  gggenes, chloroplot, grey, igv.
        gene_label_fontsize: {int} gene label fontsize.
        gene_label_color: {str} gene label color.
        show_legend: {bool} show fig legend.
        legend_size: {int} legend size.
        legend_postion: {tuple} such as (x, y), x=0-1, y=0-1.
        start: {str} initial feature, such as, ND1, ND2, ND3, ND4, ND4L, ND5, ND6,
                     COX1, COX2, COX3, ATPase6, ATPase8, Cytb, tRNA-His, tRNA-Pro,
                     tRNA-Thr, tRNA-Trp, tRNA-Met, tRNA-Asp, tRNA-Ala, tRNA-Gln,
                     tRNA-Ile, tRNA-Arg, tRNA-Tyr, tRNA-Phe, tRNA-Lys, tRNA-Gly,
                     tRNA-Asn, tRNA-Leu, tRNA-Glu, tRNA-Val, tRNA-Cys, tRNA-Ser,
                     12S rRNA, 16S rRNA, D-loop.
        height: {0-1} box height.
        species_fontsize: {int} species name fontsize.
        species_offset: {float} species name offset.
        axes: {matplotlib.projections.polar.PolarAxes}
        add_id: {bool} Species add to accession id from NCBI.
    """
    
    if not (isinstance(files, list)) and (not isinstance(files, tuple)):
        files = [files]
    
    genomes = []
    for file in reversed(files):
        features = get_features(file, abbr=abbr, isfilename2species=isfilename2species, colors=colors, start=start)
        genomes.append(features)
    
    if axes == None:
        fig, ax= plt.subplots(1, 1, figsize=(20, len(genomes)/2))
    else:
        ax = axes
    ax.set_ylim(0,len(genomes))
    
    xlimmax = 0
    for y, features in enumerate(genomes):
        for i, f in enumerate(features):
            features[i].name = FullName2AbbrName.get(f.name, f.name)
            
        if features[1].location.strand == -1:
            x = get_box_param(features[1].name)[2]+0.01
        else:
            x=0.01
            
        if add_id:
            ax.text(species_offset, y+0.5, s=features[0].name + " " + features[0].accession, style='italic', size=species_fontsize, ha="right", va="center")
        else:
            ax.text(species_offset, y+0.5, s=features[0].name, style='italic', size=species_fontsize, ha="right", va="center")
        
        for i, feature in enumerate(features[1:]):
            head_length, tail_length, box_width, text_offset= get_box_param(feature.name)
            if feature.location.strand == 1:
                arrow = FancyArrow(x, y+0.5, tail_length, 0, width=height, fc=feature.color, ec='black', head_width=height, head_length=head_length)
                ax.text(x+text_offset, y+0.5, feature.name, size=gene_label_fontsize, color=gene_label_color, ha="left", va="center",clip_on=False, clip_path=None)
                if i+1 < len(features[1:]):
                    if features[1:][i+1].location.strand == -1:
                        x += get_box_param(features[1:][i+1].name)[2] + get_box_param(features[1:][i].name)[2]
                    else:
                        x += get_box_param(features[1:][i].name)[2]
            else:
                arrow = FancyArrow(x, y+0.5, -tail_length, 0, width=height, fc=feature.color, ec='black', head_width=height, head_length=head_length)
                ax.text(x-text_offset, y+0.5, feature.name, size=gene_label_fontsize, color=gene_label_color, ha="right", va="center", clip_on=False, clip_path=None)
                if i+1 < len(features[1:]):
                    if features[1:][i+1].location.strand == -1:
                        x +=  get_box_param(features[1:][i+1].name)[2]
            ax.add_patch(arrow)
            if xlimmax < x:
                xlimmax = x
                
    ax.set_xlim(0, xlimmax+0.01)           
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    if show_legend:
        legend_elements = [Line2D([0], [0], color='black',marker='>', markersize=legend_size, 
                                  linestyle='-', linewidth=0, markerfacecolor='white', 
                                  markeredgecolor='black', label="Forward"),
                           Line2D([0], [0], color='black',marker='<', markersize=legend_size, 
                                  linestyle='-', linewidth=0, markerfacecolor='white', 
                                  markeredgecolor='black', label="Reverse")]
        if colors == None:
            legend_elements.extend(MTColors_legends['MITOFISH'])
            ncol = 12
        elif isinstance(colors, str):
            legend_elements.extend(MTColors_legends.get(colors.upper(), 'MITOFISH'))
            ncol=12
        elif isinstance(colors, dict):
            legend_elements.extend([Patch(facecolor=colors[i], edgecolor='black', label=i) for i in colors if i!="source"])
            ncol = 12
            
        ax.legend(handles=legend_elements,
                  loc='upper right', 
                  bbox_to_anchor=legend_postion, 
                  ncol=12,
                  #frameon=False,
                  #shadow=False,
                  prop={'size': legend_size},
                  title='')
                  
    if output!=None:
        fig.savefig(output, bbox_inches='tight')
    if axes != None:
        return ax
    else:
        return fig, ax

def add_tag(axs=None, 
            by_row=True, 
            fontsize=20, 
            fontstyle='bold',
            fontpostiton=(-0.15, 1.1), 
            labels = ["A", "B", "C", "D", "E", "F",
                      "G", "H", "I", "J", "K", "L",
                      "M", "N", "O", "P", "Q", "I",
                      "S", "T", "U", "V", "W", "X",
                      "Y", "Z"]):
    """
    Descripton:
        Add tags to multiple axes.
        
    Parameters：
       axs: {array} a include axes array.
       by_row: {bool}
       fontsize {float}
       fontstyle: {str} font style. such as plain, bold, italic, bold.italic
       fontpostiton: {tuple} font postition.
       labels: a set labels, such as, ["A", "B", "C", ...], ["I", "II", "III", ...].
    
    """
    fontstyles = {"plain": ("normal", "normal"),
     "bold": ("bold", "normal"),
     "italic": ("normal", "italic"),
     "bold.italic": ("bold", "italic")
    }
    
    axs_new= []
    if len(axs.shape) == 1:
        axs_new = axs
    else:
        if by_row:
            for x in range(0, axs.shape[0]):
                for y in range(0, axs.shape[1]):
                    axs_new.append(ax[(x,y)])
        else:
            for x in range(0, axs.shape[1]):
                for y in range(0, axs.shape[0]):
                    axs_new.append(ax[(y,x)])
    for n, ax in enumerate(axs_new):
        
        
        ax.text(*fontpostiton, labels[n], transform=ax.transAxes, fontsize=fontsize, 
                fontweight=fontstyles.get(fontstyle, ("bold", "normal"))[0],
                fontstyle=fontstyles.get(fontstyle, ("bold", "normal"))[1],
                verticalalignment='top')
    return ax
