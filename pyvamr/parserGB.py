#!/usr/bin/env python
# coding: utf-8

import re
import os
import sys
import logging
from Bio import SeqIO, Entrez
from .config import CommonNamesDict, MTColors
from Bio.SeqFeature import CompoundLocation, SimpleLocation

class Feature:
    def __init__(self, name, location, type, color, join=False, mtgenome=None, accession=None, file=None):
        self.name = name
        self.location = location
        self.type = type
        self.color = color
        self.join = join
        self.mtgenome = mtgenome
        self.accession = accession
        self.file = file
    def __repr__(self):
            return str(self.name)
    def __str__(self):
            return self.__repr__()

def search_name(gene_name):
    gene_name = gene_name.upper()
    status = False
    if CommonNamesDict.get(gene_name) == None:
        for p in CommonNamesDict:
            if bool(re.search(p+".", gene_name)):
                status = True
                return CommonNamesDict.get(p)
    if status == False:
        return gene_name
    else:
        return CommonNamesDict.get(gene_name)

def is_repeat(features, location):
    status = False
    for i in features:
        if i.location == location:
            status = True
    return status

def rotate_seq(seq, index):
    k = len(seq) - index
    k = k % len(seq)
    return seq[-k:] + seq[:-k]

def get_genbank_from_ncbi(accession):
    Entrez.email = "thecgs001.foxmail.com"
    gb_text = Entrez.efetch(db="Nucleotide", id=accession, rettype='gb')
    return gb_text
    
def get_species_name(string, abbr=True):
    species = re.sub("_", " ", string)
    if abbr:
        if bool(re.search(' x ', species)):
            tmp = species.split(' x ')
            species1 = tmp[0].split(' ')[0][0] + '. ' +  tmp[0].split(' ')[1]
            species2 = tmp[1].split(' ')[0][0] + '. ' +  tmp[1].split(' ')[1]
            species = species1 + ' x ' + species2
        else:
            if len(species.split(' ')) == 2:
                species = species.split(' ')[0][0] + '. ' + species.split(' ')[1]
            if len(species.split(' ')) == 3:
                species = species.split(' ')[0][0] + '. ' + species.split(' ')[1] + ' ' + ' '.join(species.split(' ')[2:])
    return species

def reinit_features(features, start = "tRNA-Phe"):
    logger = logging.getLogger(__name__) 
    logger.setLevel(logging.DEBUG)
    #logging.basicConfig(level=logging.DEBUG)
    
    features_new = [features[0]]
    
    value = None
    if isinstance(start, str):
        for feature in features[1:]:
            if feature.name == start:
                value = feature.location.start
                break
    if value == None:
        logger.warning(features[0].file+" does not contain " + start + ". Therefore, the rotation operation was not performed. Please check.")
        return features
        
    _status = True
    for index, feature in enumerate(features[1:]):
        #print(feature, feature.location, feature.name)
        if feature.name == start:
            _status = False
        if _status:
            location = SimpleLocation(start=feature.location.start + len(features[0].location)-value,
                                      end=feature.location.end + len(features[0].location)-value,
                                      strand=feature.location.strand)
            features_new.append(Feature(feature.name, location, feature.type, feature.color, feature.join))
        else:
            location = SimpleLocation(start=feature.location.start-value,
                                      end=feature.location.end-value,
                                      strand=feature.location.strand)
            features_new.append(Feature(feature.name, location, feature.type, feature.color, feature.join))
    tmp = [features_new[0]]
    tmp.extend(sorted(features_new[1:], key=lambda x:x.location.start))
    f = Feature(features[0].name, features[0].location, features[0].type, features[0].color, features[0].join, rotate_seq(features[0].mtgenome, value), accession=features[0].accession, file=features[0].file)
    features_new = [f]
    index = 0
    for feature in tmp[1:]:
        if (features_new[index].name != feature.name) or (features_new[index].join==False) or (feature.join==False):
            features_new.append(feature)
            index += 1
        else:
            location = SimpleLocation(start=features_new[index].location.start,
                                      end=feature.location.end,
                                      strand=feature.location.strand)
            features_new[index].location = location

    return features_new
    

def get_features(file, abbr=False, colors=None, isfilename2species=False, start="tRNA-Phe"):
    if colors == None:
        colors = MTColors['MITOFISH']
    elif isinstance(colors, str):
        colors = MTColors.get(colors.upper(), 'MITOFISH')
            
    features = []
    
    if os.path.exists(file):
        handle = SeqIO.parse(file, 'genbank')
    else:
        handle = SeqIO.parse(get_genbank_from_ncbi(file), 'genbank') 
    #print(handle)
    for record in handle:
        mtgenome = record.seq.upper()
        accession = record.id
        #print(file, record.id)
        for i in record.features:
            if len(list(i.qualifiers.values())) != 0:
                if "product" in i.qualifiers:
                    gene_name = i.qualifiers['product'][0]
                    gene_name =  search_name(gene_name)
                    
                    #print("gene_product", gene_name)
                    if gene_name not in CommonNamesDict and "gene" in i.qualifiers:
                        gene_name = i.qualifiers['product'][0]
                        gene_name =  search_name(gene_name)
                        
                        if gene_name not in CommonNamesDict and "note" in i.qualifiers:
                            gene_name = i.qualifiers['note'][0]
                            gene_name =  search_name(gene_name)
                        
                elif "gene" in i.qualifiers:
                    gene_name = i.qualifiers['gene'][0]
                    gene_name =  search_name(gene_name)
                    if gene_name not in CommonNamesDict and "note" in i.qualifiers:
                        gene_name = i.qualifiers['note'][0]
                        gene_name =  search_name(gene_name)
                
                elif "note" in i.qualifiers:
                    gene_name = i.qualifiers["note"][0]
                    gene_name =  search_name(gene_name)
                    
                elif "organism" in i.qualifiers:
                    gene_name = i.qualifiers["organism"][0]
                    #gene_name =  search_name(gene_name)
                    
                else:
                    #print(i.qualifiers)
                     pass
            else:
                gene_name = i.type
                gene_name =  search_name(gene_name)

            gene_name = CommonNamesDict.get(gene_name.upper(), gene_name)
            #print(gene_name)
            
            if i.type ==  "source":
                if isfilename2species:
                    species_name = os.path.splitext(os.path.basename(file))[0]
                species_name = i.qualifiers["organism"][0]
                species_name = get_species_name(species_name, abbr=abbr)
                
                features.append(Feature(name=species_name, location=i.location, type=i.type, color=colors.get('source', 'gray'),
                                        mtgenome=mtgenome, accession=accession, file=file))
                                        
            elif i.type in ['rRNA', 'tRNA', 'D_loop', 'D-loop']:
                if gene_name in ['tRNA-His', 'tRNA-Pro', 'tRNA-Thr', 'tRNA-Trp', 'tRNA-Met', 'tRNA-Asp', 'tRNA-Ala', 'tRNA-Gln',
                                 'tRNA-Ile', 'tRNA-Arg', 'tRNA-Tyr', 'tRNA-Phe', 'tRNA-Lys', 'tRNA-Gly', 'tRNA-Asn', 'tRNA-Leu',
                                 'tRNA-Glu', 'tRNA-Val', 'tRNA-Cys', 'tRNA-Ser', '12S rRNA', '16S rRNA', "D-loop"]:
                    if isinstance(i.location, CompoundLocation):
                        for location in i.location.parts:
                            if not is_repeat(features=features, location=location):
                                features.append(Feature(name=gene_name, location=location, type=i.type, color=colors.get(gene_name, 'gray'),join=True))
                    else:
                        if not is_repeat(features=features, location=i.location):
                            features.append(Feature(name=gene_name, location=i.location, type=i.type, color=colors.get(gene_name, 'gray')))
                        
            elif i.type in ['CDS', 'gene']:
                if gene_name in ['ND1', 'ND2', 'ND3', 'ND4L', 'ND4', 'ND5', 'ND6', 'COX1', 'COX2', 'COX3', 'ATPase6', 'ATPase8', 'Cytb']:
                    if isinstance(i.location, CompoundLocation):
                        for location in i.location.parts:
                            if not is_repeat(features=features, location=location):
                                features.append(Feature(name=gene_name, location=location, type=i.type, color=colors.get(gene_name, 'gray'),join=True))
                    else:
                        if not is_repeat(features=features, location=i.location):
                            features.append(Feature(name=gene_name, location=i.location, type=i.type, color=colors.get(gene_name, 'gray')))
                        
            elif i.type == 'misc_feature':
                if gene_name in ['tRNA-His', 'tRNA-Pro', 'tRNA-Thr', 'tRNA-Trp', 'tRNA-Met', 'tRNA-Asp', 'tRNA-Ala', 'tRNA-Gln',
                                 'tRNA-Ile', 'tRNA-Arg', 'tRNA-Tyr', 'tRNA-Phe', 'tRNA-Lys', 'tRNA-Gly', 'tRNA-Asn', 'tRNA-Leu',
                                 'tRNA-Glu', 'tRNA-Val', 'tRNA-Cys', 'tRNA-Ser', '12S rRNA', '16S rRNA', "D-loop"]:
                    if isinstance(i.location, CompoundLocation):
                        for location in i.location.parts:
                            features.append(Feature(name=gene_name, location=location, type=i.type, color=colors.get(gene_name, 'gray'),join=True))
                    else:
                        features.append(Feature(name=gene_name, location=i.location, type=i.type, color=colors.get(gene_name, 'gray')))
            else:
                pass
    
    #print(features)
    if features == []:
        logger = logging.getLogger(__name__) 
        logger.setLevel(logging.DEBUG)
        logger.warning(file + " Gene ID error!")
        sys.exit()
    
    res = [features[0]]
    
    res.extend(sorted(features[1:], key=lambda x:x.location.start))
    
    if start !=None:
        res = reinit_features(res, start = start)
    return res
