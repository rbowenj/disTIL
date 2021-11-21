from numpy.core.numeric import NaN
import pandas as pd
from argparse import ArgumentParser
import re


def create_arg_parser():
    # Create argument parser
    description = "Filter neoepitopes predicted by pVACfuse."
    argparser = ArgumentParser(description = description)
    argparser.add_argument('pvacfuse_report', type=str,
                            help='Aggregated pVACfuse epitope report.')
    argparser.add_argument('star_fusions', type=str, 
                            help='STAR-Fusion TSV of fusions for this patient.')
    argparser.add_argument('-n_junction_reads', type=int, 
                            help='Junction reads count filter threshold. Neoepitopes derived from fusions with junction read counts below this threshold iwll be discarded.')
    argparser.add_argument('-n_spanning_frags', type=int, 
                            help='Spanning fragment count filter threshold. Neoepitopes derived from fusions with spanning fragment counts below this threshold iwll be discarded.')
    # argparser.add_argument('patient_id', type=str,
    #                         help='Patient ID.')
    
    return argparser

def read_file(filename):
    df = pd.read_csv(filename, sep='\t', header=0)
    return df

def reads_filter(df, n, colname):
    # Only keep cols with at least one HLA match (class 1 and/or 2)
    return df[df[colname] >= n]

if __name__ == '__main__':
    argparser = create_arg_parser()
    args = argparser.parse_args()

    # pvacfuse_df = read_file(args.pvacfuse_report)
    # star_df = read_file(args.star_fusions)
    pvacfuse_df = pd.read_csv(args.pvacfuse_report, sep="\t", header=0)
    star_df = pd.read_csv(args.star_fusions, sep="\t", header=0)

    # star_df_sub = star_df[['#FusionName', 'JunctionReadCount', 'SpanningFragCount']]

    star_df['left_gene'] = None
    star_df['right_gene'] = None
    star_df['left_breakpoint'] = None
    star_df['right_breakpoint'] = None

    pvacfuse_df['left_gene'] = None
    pvacfuse_df['right_gene'] = None
    pvacfuse_df['breakpoints'] = None
    pvacfuse_df['left_breakpoint'] = None
    pvacfuse_df['right_breakpoint'] = None

    for i, r in star_df.iterrows():
        gene_split = r['#FusionName'].split('--')
        new_gene = '-'.join(gene_split)

        star_df.at[i,'#FusionName'] = new_gene
        star_df.at[i, 'left_gene'] = gene_split[0]
        star_df.at[i, 'right_gene'] = gene_split[1]

        left = re.search(r'\:(\d+)', r['LeftBreakpoint'])
        star_df.at[i, 'left_breakpoint'] = left.group(1)

        right = re.search(r'\:(\d+)', r['RightBreakpoint'])
        star_df.at[i, 'right_breakpoint'] = right.group(1)

    for i, neoep in pvacfuse_df.iterrows():
        gene_split = neoep['Gene Name'].split('-')
        pvacfuse_df.at[i, 'left_gene'] = gene_split[0]
        pvacfuse_df.at[i, 'right_gene'] = gene_split[1]
        pvacfuse_df.at[i, 'breakpoints'] = [x.strip() for x in str(neoep['Stop']).split('/') + str(neoep['Start']).split('/')]
        # pvacfuse_df.at[i, 'left_breakpoint'] = str(neoep['Stop']).split('/')[0].strip()
        # pvacfuse_df.at[i, 'right_breakpoint'] = str(neoep['Start']).split('/')[1].strip()

        for j, fusion in star_df.iterrows():
            if fusion['left_breakpoint'] in pvacfuse_df.iloc[i]['breakpoints'] and fusion['right_breakpoint'] in pvacfuse_df.iloc[i]['breakpoints']:
                pvacfuse_df.at[i, 'left_breakpoint'] = fusion['left_breakpoint']
                pvacfuse_df.at[i, 'right_breakpoint'] = fusion['right_breakpoint']
    
    star_df = star_df[['left_gene', 'right_gene', 'left_breakpoint', 'right_breakpoint', 'JunctionReadCount', 'SpanningFragCount']]
    merged_df = pvacfuse_df.merge(star_df, on=['left_gene', 'right_gene', 'left_breakpoint', 'right_breakpoint'], how='left')
    
    if args.n_junction_reads:
        merged_df = filter(merged_df, args.n_junction_reads, 'JunctionReadCount')

    if args.n_spanning_frags:
        merged_df = filter(merged_df, args.n_spanning_frags, 'SpanningFragCount')

    # Subset columns
    merged_df = merged_df.drop(['left_gene', 'right_gene', 'left_breakpoint', 'right_breakpoint', 'breakpoints'], axis=1)
    merged_df = merged_df.rename(columns = {'JunctionReadCount': '# Junction reads', 'SpanningFragCount': '# Spanning fragments'})

    infile_base = args.pvacfuse_report.split('/')[-1].split('.tsv')[0]
    outfile = infile_base + '.junctionFiltered.tsv'
    merged_df.to_csv(outfile, sep='\t', header=True)