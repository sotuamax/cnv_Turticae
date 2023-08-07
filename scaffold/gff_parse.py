
import pandas as pd
import argparse
import statistics
import collections
import os
import subprocess

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="Parse gff and clean gene coordinates on contigs for its chromosome assignment.")
    parser.add_argument("target_gff", help = "gff file for target genome. ")
    parser.add_argument("ref_gff", help = "gff file for reference genome. ")
    parser.add_argument("-vis", "--vis", required = False, help = "path for the visualization tool. ")
    parser.add_argument("-O", "--output", help="output prefix")
    args = parser.parse_args()
    return args

def gff_loc(args):
    '''parse reference gff file and retrieve gene locations. '''
    ref = args.ref_gff 
    refgff = pd.read_table(ref, sep = "\t", header = None)
    refgff = refgff.rename(columns = {0:"chromosome", 1:"source", 2:"feature", 3:"start", 4:"end", 5:"quality", 6:"strand", 7:".", 8:"attribute"})
    refgene = refgff[refgff["feature"] == "gene"].copy()
    refgene["gene"] = refgene.attribute.str.split(";", expand = True)[0].str.split("gene:", expand = True)[1]
    refgene == refgene[["gene", "chromosome", "start", "end", "strand"]]
    return refgene

def parse_gff(args):
    '''parse gff file to only contain consecutive genes.'''
    gff = args.gff
    loc = gff_loc(args)
    output = args.output
    ###
    gffmap = pd.read_table(gff, sep = "\t", header = None)
    gffmap = gffmap.rename(columns = {0:"contig", 1:"source", 2:"feature", 3:"tig_start", 4:"tig_end", 5:"quality", 6:"tig_strand", 7:".", 8:"attribute"})
    gffgene = gffmap[gffmap["feature"] == "gene"].copy() # only gene rows
    gffgene["gene"] = gffgene.attribute.str.split(";", expand = True)[0].str.split("gene:", expand = True)[1]
    ###
    gffloc = gffgene.merge(loc, on = "gene", how = "left")
    gffloc = gffloc.sort_values(by = ["chromosome", "start"])
    gffloc["gorder"] = range(1, len(gffloc)+1)
    ###
    gffmap = consec_gtf(gffloc)
    gffmap = gffmap.sort_values(by = ["contig", "tig_start"], ascending = True)
    gffmap["torder"] = range(1, len(gffmap)+1)
    gffmap_chr = gffmap[gffmap.chromosome.str.contains("chromosome")]
    gffmap_chr.to_csv(output+".gff", sep = "\t", index = False)
    return gffmap

def get_freq(args, gffmap):
    """Genes frequency on contigs matched to chromosome. """
    ###
    output = args.output
    tigs = sorted(set(gffmap["contig"]))
    tmap_dict = {t:collections.Counter(gffmap[gffmap["contig"] == t]["chromosome"]) for t in tigs}
    fw = open(output+".freq", "w")
    fw.write("contig\tchromosome\tfreq\n")
    for t in tmap_dict:
        for c in tmap_dict[t]:
            fw.write("{}\t{}\t{}\n".format(t, c, tmap_dict[t][c]))
    fw.close()
    ### frequency for all reference genome including chromosome and scaffold
    freq_df = pd.read_table(output+".freq", sep = "\t", header = 0)
    df_container = list()
    for t in tmap_dict:
        freq_sub = freq_df[freq_df["contig"] == t].copy()
        freqsum = sum(freq_sub["freq"])
        freq_sub["prop"] = round(freq_sub["freq"]/freqsum*100, 2)
        df_container.append(freq_sub)
    ###
    prop_df_raw = pd.concat(df_container, axis = 0)
    prop_df_raw.to_csv(output + ".freq", sep = "\t", index = False)
    ### frequency for reference genome only with chromosome
    freq_df = pd.read_table(output + ".freq", sep = "\t", header = 0)
    df_container = list()
    for t in tmap_dict:
        freq_sub = freq_df[(freq_df["contig"] == t) & (freq_df["chromosome"].str.contains("chromosome"))].copy()
        freqsum = sum(freq_sub["freq"])
        freq_sub["prop"] = round(freq_sub["freq"]/freqsum*100, 2)
        df_container.append(freq_sub)
    ###
    prop_df_chr = pd.concat(df_container, axis = 0)
    #prop_df_chr.to_csv(output + ".freq", sep = "\t", index = False)
    ###
    return prop_df_chr

def consec_gtf(gtf_qual):
    """if gff file with genes in consecutive order. """
    #gtf_qual = gtf_qual[~gtf_qual["chromosome"].isna()]
    tigs = sorted(set(gtf_qual["contig"]))
    cons_tig = list()
    for c in tigs:
        tig_gorder = list(gtf_qual[gtf_qual["contig"]==c]["gorder"].values[:])
        cons_num = consecutive_start(tig_gorder)
        subgtf = gtf_qual[gtf_qual["gorder"].isin(cons_num)]
        cons_tig.append(subgtf)
    tig_df = pd.concat(cons_tig, axis = 0)
    return tig_df
    #tig_df.to_csv(args.output + ".txt", sep = "\t", index = False)

def consecutive_start(l):
    """Function for consecutive order (less than 5 difference)"""
    sl = sorted(l)
    ll = len(l)
    sl_list = list()
    for i in range(ll-1):
        if sl[i+1] - sl[i] < 5:
            sl_list.append(sl[i+1])
            sl_list.append(sl[i])
    sl_list = sorted(set(sl_list))
    return sl_list

def assign_chr(args, freq_df, gffmap):
    '''contig and its potential mapped chromosome'''
    output = args.output
    ##
    tigs = sorted(set(freq_df["contig"]))
    fw = open(output+".asn", "w")
    fw.write("contig\tchromosome\tpgene\n") # pgene: potentially contained genes in the contig
    for t in tigs: # contigs have at least 3 genes on it for consideration
        freq_sub = freq_df[(freq_df["contig"] == t) & (freq_df["prop"] > 80) & (freq_df["freq"] > 1)]
        if len(freq_sub) > 0:
            pchr = freq_sub["chromosome"].values[0]
            pfreq = freq_sub["freq"].values[0]
            fw.write("{}\t{}\t{}\n".format(t, pchr, pfreq))
    fw.close()
    ###
    df_asn = pd.read_table(output+".asn", sep = "\t", header = 0)
    strand_list = list()
    order_list = list()
    for row in df_asn.itertuples():
        t = row.contig
        m = row.chromosome
        submap = gffmap[(gffmap["contig"] == t) & (gffmap["chromosome"] == m)].copy()
        suborder = statistics.median(submap["gorder"])
        order_list.append(suborder)
        submap['strand_measure'] = (submap["tig_strand"] == submap["strand"]).astype(int)
        col_dict = collections.Counter(submap["strand_measure"])
        if col_dict[1] > col_dict[0]:
            strand_list.append("+")
        else:
            strand_list.append("-")
    df_asn["strand"] = strand_list
    df_asn["tig_median"] = order_list
    df_asn = df_asn.sort_values(by = ["chromosome", "tig_median"])
    df_asn.to_csv(output+".asn", sep = "\t", index = False)
    return df_asn

def vis_gff(args):
    """visualize parsed gff file with gene order on chromosome and target contigs. """
    output = args.output
    vis = args.vis_path
    command = f"Rscript {vis} -gff {output}.gff -O {output}"
    subprocess.call(command, shell = True)

def main():
    """Run all functions. """
    args = args_parser()
    print("assign source-chromosome for contigs ......")
    parsed_gff = parse_gff(args)
    gff_freq_chr = get_freq(args, parsed_gff)
    tig_map = assign_chr(args, gff_freq_chr, parsed_gff)
    if args.vis != None: 
        vis_gff(args)

##############
### Run it ###
##############

if __name__ == "__main__":
    main()
